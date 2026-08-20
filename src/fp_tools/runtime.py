"""Provision and select pinned external-tool runtimes for fp-tools.

Normal Python and desktop installations use the managed runtime automatically
when a workflow needs command-line genomics programs.  System programs and the
complete container remain explicit advanced backends.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import platform
import shutil
import subprocess
import sys
import tarfile
import tempfile
import time
import urllib.parse
import urllib.request
from contextlib import contextmanager
from dataclasses import dataclass
from importlib import resources
from pathlib import Path
from typing import Iterable, Sequence

from fp_tools import __version__


RUNTIME_MODES = ("auto", "managed", "system", "container")
OCI_ACCEPT = (
    "application/vnd.oci.image.manifest.v1+json,"
    "application/vnd.docker.distribution.manifest.v2+json"
)


class RuntimeProvisionError(RuntimeError):
    """Raised when the managed external-tool runtime cannot be prepared."""


@dataclass(frozen=True)
class RuntimeActivation:
    mode: str
    component: str
    prefix: Path | None = None
    distro: str | None = None


def add_runtime_argument(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "--runtime",
        choices=RUNTIME_MODES,
        default=os.environ.get("FP_TOOLS_RUNTIME", "auto"),
        help=(
            "External-tool runtime: auto/managed provisions the pinned fp-tools "
            "runtime, system uses PATH, and container uses the complete image "
            "(default: auto)."
        ),
    )


def platform_key() -> str:
    system = platform.system().lower()
    machine = platform.machine().lower()
    arch = {
        "amd64": "x86_64",
        "x64": "x86_64",
        "x86_64": "x86_64",
        "aarch64": "arm64",
        "arm64": "arm64",
    }.get(machine, machine)
    if system == "darwin":
        return f"macos-{arch}"
    if system == "windows":
        return f"windows-{arch}"
    if system == "linux":
        return f"linux-{arch}"
    raise RuntimeProvisionError(f"Unsupported runtime platform: {system}-{machine}")


def runtime_cache_root() -> Path:
    override = os.environ.get("FP_TOOLS_RUNTIME_CACHE")
    if override:
        return Path(override).expanduser().resolve()
    if os.name == "nt":
        base = Path(os.environ.get("LOCALAPPDATA", Path.home() / "AppData" / "Local"))
        return base / "fp-tools" / "runtimes"
    if platform.system() == "Darwin":
        return Path.home() / "Library" / "Caches" / "fp-tools" / "runtimes"
    base = Path(os.environ.get("XDG_CACHE_HOME", Path.home() / ".cache"))
    return base / "fp-tools" / "runtimes"


def load_runtime_manifest() -> dict:
    manifest = resources.files("fp_tools.resources").joinpath("runtime_manifest.json")
    return json.loads(manifest.read_text(encoding="utf-8"))


def _runtime_mode(value: str | None) -> str:
    mode = str(value or os.environ.get("FP_TOOLS_RUNTIME", "auto")).lower()
    if mode not in RUNTIME_MODES:
        raise RuntimeProvisionError(
            f"Unknown runtime mode {mode!r}; choose {', '.join(RUNTIME_MODES)}"
        )
    return "managed" if mode == "auto" else mode


def _runtime_spec(component: str, target_platform: str | None = None) -> dict:
    manifest = load_runtime_manifest()
    target = target_platform or platform_key()
    try:
        artifact = manifest["artifacts"][target][component]
    except KeyError as exc:
        raise RuntimeProvisionError(
            f"Managed runtime component {component!r} is not available for {target}."
        ) from exc
    spec = {
        "schema": manifest["schema"],
        "runtime_version": manifest["runtime_version"],
        "commands": manifest["components"][component]["commands"],
        **artifact,
    }
    for key in ("repository", "release_base_url"):
        if key in manifest:
            spec[key] = manifest[key]
    return spec


def _request_json(url: str, headers: dict[str, str] | None = None) -> tuple[dict, dict]:
    request = urllib.request.Request(url, headers=headers or {})
    with urllib.request.urlopen(request, timeout=60) as response:
        return json.load(response), dict(response.headers)


def _request_text(url: str) -> str:
    request = urllib.request.Request(url)
    with urllib.request.urlopen(request, timeout=60) as response:
        return response.read().decode("utf-8")


def _oci_bearer_token(repository: str) -> str:
    registry, image = repository.split("/", 1)
    query = urllib.parse.urlencode(
        {"service": registry, "scope": f"repository:{image}:pull"}
    )
    payload, _ = _request_json(f"https://{registry}/token?{query}")
    token = payload.get("token") or payload.get("access_token")
    if not token:
        raise RuntimeProvisionError(f"No anonymous pull token returned for {repository}")
    return str(token)


def _resolve_oci_layer(repository: str, tag: str) -> tuple[str, int, str]:
    registry, image = repository.split("/", 1)
    token = _oci_bearer_token(repository)
    headers = {"Authorization": f"Bearer {token}", "Accept": OCI_ACCEPT}
    manifest, _ = _request_json(
        f"https://{registry}/v2/{image}/manifests/{urllib.parse.quote(tag, safe='')}",
        headers,
    )
    layers = manifest.get("layers") or []
    if len(layers) != 1:
        raise RuntimeProvisionError(
            f"Runtime artifact {repository}:{tag} must contain exactly one layer."
        )
    layer = layers[0]
    digest = str(layer.get("digest", ""))
    if not digest.startswith("sha256:"):
        raise RuntimeProvisionError(f"Runtime artifact {repository}:{tag} has no SHA-256 digest.")
    size = int(layer.get("size", 0))
    url = f"https://{registry}/v2/{image}/blobs/{digest}"
    return url, size, digest.split(":", 1)[1]


def _resolve_runtime_artifact(spec: dict) -> tuple[str, int, str]:
    filename = spec.get("filename")
    release_base_url = spec.get("release_base_url")
    if filename and release_base_url:
        url = f"{str(release_base_url).rstrip('/')}/{urllib.parse.quote(str(filename))}"
        checksum = _request_text(url + ".sha256").strip().split()[0].lower()
        if len(checksum) != 64 or any(char not in "0123456789abcdef" for char in checksum):
            raise RuntimeProvisionError(f"Runtime artifact {filename} has no valid SHA-256 checksum.")
        return url, 0, checksum
    return _resolve_oci_layer(spec["repository"], spec["tag"])


def _download(url: str, destination: Path, expected_size: int, expected_sha256: str) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    partial = destination.with_suffix(destination.suffix + ".part")
    existing = partial.stat().st_size if partial.exists() else 0
    headers = {"Range": f"bytes={existing}-"} if existing else {}
    if "ghcr.io/" in url:
        repository = load_runtime_manifest()["repository"]
        headers["Authorization"] = f"Bearer {_oci_bearer_token(repository)}"
    request = urllib.request.Request(url, headers=headers)
    try:
        with urllib.request.urlopen(request, timeout=120) as response:
            resumed = bool(existing and getattr(response, "status", None) == 206)
            if existing and not resumed:
                existing = 0
            mode = "ab" if resumed else "wb"
            handle = partial.open(mode)
            try:
                total = expected_size or int(response.headers.get("Content-Length", 0)) + existing
                downloaded = existing
                last_report = -1
                while True:
                    chunk = response.read(1024 * 1024)
                    if not chunk:
                        break
                    handle.write(chunk)
                    downloaded += len(chunk)
                    percent = int(downloaded * 100 / total) if total else 0
                    if percent >= last_report + 5:
                        print(f"Preparing fp-tools runtime: {percent}%", flush=True)
                        last_report = percent
            finally:
                handle.close()
    except Exception as exc:
        raise RuntimeProvisionError(
            "Runtime download was interrupted; rerun the command to resume."
        ) from exc
    if expected_size and partial.stat().st_size != expected_size:
        raise RuntimeProvisionError(
            f"Runtime download size mismatch: expected {expected_size}, got {partial.stat().st_size}."
        )
    digest = hashlib.sha256()
    with partial.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    if digest.hexdigest() != expected_sha256:
        partial.unlink(missing_ok=True)
        raise RuntimeProvisionError("Runtime checksum verification failed; the partial download was removed.")
    os.replace(partial, destination)


def _safe_extract(archive: Path, destination: Path) -> None:
    destination_resolved = destination.resolve()
    with tarfile.open(archive, "r:gz") as handle:
        for member in handle.getmembers():
            target = (destination / member.name).resolve()
            if destination_resolved != target and destination_resolved not in target.parents:
                raise RuntimeProvisionError(f"Unsafe path in runtime archive: {member.name}")
            if member.issym() or member.islnk():
                link_target = (target.parent / member.linkname).resolve()
                if destination_resolved != link_target and destination_resolved not in link_target.parents:
                    raise RuntimeProvisionError(f"Unsafe link in runtime archive: {member.name}")
        try:
            handle.extractall(destination, filter="fully_trusted")
        except TypeError:  # Python versions before tarfile extraction filters.
            handle.extractall(destination)


@contextmanager
def _install_lock(path: Path, timeout: float = 600.0):
    path.parent.mkdir(parents=True, exist_ok=True)
    deadline = time.monotonic() + timeout
    while True:
        try:
            fd = os.open(path, os.O_CREAT | os.O_EXCL | os.O_WRONLY)
            os.write(fd, f"{os.getpid()}\n".encode())
            os.close(fd)
            break
        except FileExistsError:
            try:
                if time.time() - path.stat().st_mtime > 7200:
                    path.unlink()
                    continue
            except FileNotFoundError:
                continue
            if time.monotonic() >= deadline:
                raise RuntimeProvisionError("Timed out waiting for another runtime installation.")
            time.sleep(0.25)
    try:
        yield
    finally:
        path.unlink(missing_ok=True)


def _free_space_check(path: Path, required_bytes: int) -> None:
    path.mkdir(parents=True, exist_ok=True)
    free = shutil.disk_usage(path).free
    required = max(required_bytes * 3, 512 * 1024 * 1024)
    if free < required:
        raise RuntimeProvisionError(
            f"Insufficient disk space for the fp-tools runtime: need approximately "
            f"{required / 1024**3:.1f} GiB, available {free / 1024**3:.1f} GiB."
        )


def ensure_native_runtime(component: str = "core") -> RuntimeActivation:
    target = platform_key()
    if target.startswith("windows-"):
        return ensure_wsl_runtime(component)
    spec = _runtime_spec(component, target)
    root = runtime_cache_root() / spec["runtime_version"] / target
    prefix = root / component
    ready = prefix / ".fp-tools-runtime.json"
    if ready.is_file():
        return RuntimeActivation("managed", component, prefix=prefix)
    with _install_lock(root / f".{component}.lock"):
        if ready.is_file():
            return RuntimeActivation("managed", component, prefix=prefix)
        url, size, digest = _resolve_runtime_artifact(spec)
        _free_space_check(root, size)
        archive = root / f"{component}.tar.gz"
        _download(url, archive, size, digest)
        temporary = Path(tempfile.mkdtemp(prefix=f".{component}-", dir=root))
        try:
            _safe_extract(archive, temporary)
            python = temporary / "bin" / "python"
            unpack = temporary / "bin" / "conda-unpack"
            if python.is_file() and unpack.is_file():
                result = subprocess.run([str(python), str(unpack)], check=False)
                if result.returncode:
                    raise RuntimeProvisionError("The downloaded runtime could not be relocated.")
            missing = [name for name in spec["commands"] if not (temporary / "bin" / name).exists()]
            if missing:
                raise RuntimeProvisionError(
                    f"The downloaded runtime is missing: {', '.join(missing)}"
                )
            (temporary / ".fp-tools-runtime.json").write_text(
                json.dumps(
                    {
                        "runtime_version": spec["runtime_version"],
                        "component": component,
                        "platform": target,
                        "sha256": digest,
                    },
                    indent=2,
                ),
                encoding="utf-8",
            )
            if prefix.exists():
                shutil.rmtree(prefix)
            os.replace(temporary, prefix)
        finally:
            if temporary.exists():
                shutil.rmtree(temporary, ignore_errors=True)
            archive.unlink(missing_ok=True)
    return RuntimeActivation("managed", component, prefix=prefix)


def activate_runtime(component: str = "core", mode: str | None = None) -> RuntimeActivation:
    resolved = _runtime_mode(mode)
    if resolved == "system":
        return RuntimeActivation("system", component)
    if resolved == "container":
        return RuntimeActivation("container", component)
    activation = ensure_native_runtime(component)
    if activation.prefix is not None:
        bin_dir = activation.prefix / ("Scripts" if os.name == "nt" else "bin")
        os.environ["PATH"] = str(bin_dir) + os.pathsep + os.environ.get("PATH", "")
        os.environ["FP_TOOLS_RUNTIME_PREFIX"] = str(activation.prefix)
    return activation


def _wsl_name(version: str) -> str:
    safe = "".join(ch if ch.isalnum() else "-" for ch in version).strip("-")
    return f"fp-tools-{safe}"


def _wsl_distributions() -> set[str]:
    result = subprocess.run(
        ["wsl.exe", "--list", "--quiet"], capture_output=True, text=True, check=False
    )
    text = (result.stdout or "").replace("\x00", "")
    return {line.strip() for line in text.splitlines() if line.strip()}


def _enable_wsl() -> None:
    if shutil.which("wsl.exe"):
        status = subprocess.run(["wsl.exe", "--status"], capture_output=True, check=False)
        if status.returncode == 0:
            return
    powershell = shutil.which("powershell.exe") or shutil.which("powershell")
    if not powershell:
        raise RuntimeProvisionError("Windows Subsystem for Linux is unavailable on this system.")
    command = (
        "Start-Process wsl.exe -Verb RunAs -Wait "
        "-ArgumentList '--install','--no-distribution'"
    )
    result = subprocess.run([powershell, "-NoProfile", "-Command", command], check=False)
    if result.returncode:
        raise RuntimeProvisionError("Windows did not enable WSL2.")
    raise RuntimeProvisionError(
        "WSL2 was enabled. Restart Windows if requested, then open fp-tools again; setup will resume automatically."
    )


def ensure_wsl_runtime(component: str = "core") -> RuntimeActivation:
    if os.name != "nt":
        raise RuntimeProvisionError("The WSL runtime can only be provisioned from Windows.")
    _enable_wsl()
    target = platform_key()
    spec = _runtime_spec(component, target)
    distro = _wsl_name(spec["runtime_version"])
    if distro in _wsl_distributions():
        return RuntimeActivation("managed", component, distro=distro)
    root = runtime_cache_root() / spec["runtime_version"] / target
    install_dir = root / "wsl"
    with _install_lock(root / ".wsl.lock"):
        if distro in _wsl_distributions():
            return RuntimeActivation("managed", component, distro=distro)
        url, size, digest = _resolve_runtime_artifact(spec)
        _free_space_check(root, size)
        archive = root / "fp-tools-wsl-rootfs.tar.gz"
        _download(url, archive, size, digest)
        install_dir.mkdir(parents=True, exist_ok=True)
        result = subprocess.run(
            [
                "wsl.exe",
                "--import",
                distro,
                str(install_dir),
                str(archive),
                "--version",
                "2",
            ],
            check=False,
        )
        if result.returncode:
            raise RuntimeProvisionError("Windows could not import the fp-tools WSL runtime.")
        archive.unlink(missing_ok=True)
    return RuntimeActivation("managed", component, distro=distro)


def _looks_like_url(value: str) -> bool:
    return urllib.parse.urlparse(value).scheme in {"http", "https", "ftp", "s3"}


def _replace_runtime_option(arguments: list[str], value: str) -> list[str]:
    output: list[str] = []
    skip = False
    for index, argument in enumerate(arguments):
        if skip:
            skip = False
            continue
        if argument == "--runtime":
            skip = True
            continue
        if argument.startswith("--runtime="):
            continue
        output.append(argument)
    output.extend(["--runtime", value])
    return output


def _wsl_path(distro: str, value: str) -> str:
    if _looks_like_url(value):
        return value
    path = Path(value).expanduser()
    if not path.is_absolute():
        path = (Path.cwd() / path).resolve()
    result = subprocess.run(
        ["wsl.exe", "--distribution", distro, "--exec", "wslpath", "-a", str(path)],
        capture_output=True,
        text=True,
        check=False,
    )
    if result.returncode:
        raise RuntimeProvisionError(f"Could not translate Windows path for WSL: {value}")
    return result.stdout.strip()


def translate_flag_paths(
    arguments: Sequence[str], path_flags: Iterable[str], translator
) -> list[str]:
    path_flags = set(path_flags)
    translated: list[str] = []
    index = 0
    arguments = list(map(str, arguments))
    while index < len(arguments):
        argument = arguments[index]
        if argument in path_flags and index + 1 < len(arguments):
            translated.append(argument)
            index += 1
            while index < len(arguments) and not arguments[index].startswith("--"):
                translated.append(translator(arguments[index]))
                index += 1
            continue
        matched = next((flag for flag in path_flags if argument.startswith(flag + "=")), None)
        if matched:
            translated.append(matched + "=" + translator(argument.split("=", 1)[1]))
        else:
            translated.append(argument)
        index += 1
    return translated


def run_managed_wsl_command(
    command: str,
    arguments: Sequence[str],
    path_flags: Iterable[str],
    component: str = "core",
) -> int:
    activation = ensure_wsl_runtime(component)
    assert activation.distro is not None
    translated = translate_flag_paths(
        _replace_runtime_option(list(arguments), "system"),
        path_flags,
        lambda value: _wsl_path(activation.distro or "", value),
    )
    cwd = _wsl_path(activation.distro, str(Path.cwd()))
    executable = f"/opt/conda/bin/{command}"
    result = subprocess.run(
        [
            "wsl.exe",
            "--distribution",
            activation.distro,
            "--cd",
            cwd,
            "--exec",
            executable,
            *translated,
        ],
        check=False,
    )
    return int(result.returncode)


def run_container_command(
    command: str,
    arguments: Sequence[str],
    path_flags: Iterable[str],
) -> int:
    """Run one fp-tools command in the complete public container."""

    docker = shutil.which("docker")
    if not docker:
        raise RuntimeProvisionError(
            "Container runtime requested, but Docker is not installed or not available on PATH."
        )
    cwd = Path.cwd().resolve()
    mounts: dict[Path, str] = {cwd: "/work"}

    def translate(value: str) -> str:
        if _looks_like_url(value) or value in {"hg38", "mm10"}:
            return value
        path = Path(value).expanduser()
        if not path.is_absolute():
            path = (cwd / path).resolve()
        try:
            return str(Path("/work") / path.relative_to(cwd))
        except ValueError:
            parent = path if path.exists() and path.is_dir() else path.parent
            parent = parent.resolve()
            target = mounts.setdefault(parent, f"/fp-tools-mount/{len(mounts)}")
            return target if path == parent else str(Path(target) / path.name)

    translated = translate_flag_paths(
        _replace_runtime_option(list(arguments), "system"), path_flags, translate
    )
    image = os.environ.get(
        "FP_TOOLS_CONTAINER_IMAGE", f"fp-tools:v{__version__}"
    )
    invocation = [docker, "run", "--rm"]
    if os.name != "nt" and hasattr(os, "getuid"):
        invocation.extend(["--user", f"{os.getuid()}:{os.getgid()}", "-e", "HOME=/tmp"])
    for source, target in mounts.items():
        invocation.extend(["-v", f"{source}:{target}"])
    invocation.extend(["-w", "/work", image, command, *translated])
    return int(subprocess.run(invocation, check=False).returncode)


def prepare_command_runtime(
    command: str,
    arguments: Sequence[str],
    mode: str | None,
    component: str,
    path_flags: Iterable[str],
) -> int | None:
    """Prepare a runtime or delegate the complete command when required."""

    resolved = _runtime_mode(mode)
    if resolved == "container":
        return run_container_command(command, arguments, path_flags)
    if resolved == "managed" and os.name == "nt":
        return run_managed_wsl_command(command, arguments, path_flags, component)
    activate_runtime(component, resolved)
    return None


def runtime_status() -> list[dict[str, str]]:
    manifest = load_runtime_manifest()
    target = platform_key()
    rows = []
    for component in manifest["components"]:
        available = component in manifest.get("artifacts", {}).get(target, {})
        if target.startswith("windows-"):
            installed = _wsl_name(manifest["runtime_version"]) in _wsl_distributions() if shutil.which("wsl.exe") else False
            location = _wsl_name(manifest["runtime_version"])
        else:
            location_path = runtime_cache_root() / manifest["runtime_version"] / target / component
            installed = (location_path / ".fp-tools-runtime.json").is_file()
            location = str(location_path)
        rows.append(
            {
                "component": component,
                "platform": target,
                "available": "yes" if available else "no",
                "installed": "yes" if installed else "no",
                "location": location,
            }
        )
    return rows
