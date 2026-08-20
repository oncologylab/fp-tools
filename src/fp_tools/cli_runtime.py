"""Inspect, install, or repair the managed fp-tools runtime."""

from __future__ import annotations

import argparse
import subprocess
import shutil

from fp_tools.runtime import (
    RuntimeProvisionError,
    ensure_native_runtime,
    load_runtime_manifest,
    platform_key,
    runtime_cache_root,
    runtime_status,
)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="fp-tools-runtime", description=__doc__)
    subparsers = parser.add_subparsers(dest="action", required=True)
    subparsers.add_parser("status", help="Report managed runtime availability and installation state.")
    for action in ("install", "repair"):
        child = subparsers.add_parser(action, help=f"{action.title()} a runtime component.")
        child.add_argument(
            "component",
            choices=tuple(load_runtime_manifest()["components"]),
            nargs="?",
            default="core",
        )
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    if args.action == "status":
        print("component\tplatform\tavailable\tinstalled\tlocation")
        for row in runtime_status():
            print("\t".join(row[key] for key in ("component", "platform", "available", "installed", "location")))
        return 0
    if args.action == "repair" and not platform_key().startswith("windows-"):
        target = runtime_cache_root() / load_runtime_manifest()["runtime_version"] / platform_key() / args.component
        if target.exists():
            shutil.rmtree(target)
    elif args.action == "repair":
        manifest = load_runtime_manifest()
        safe_version = "".join(
            char if char.isalnum() else "-" for char in manifest["runtime_version"]
        ).strip("-")
        if shutil.which("wsl.exe"):
            subprocess.run(
                ["wsl.exe", "--unregister", f"fp-tools-{safe_version}"],
                check=False,
            )
    try:
        activation = ensure_native_runtime(args.component)
    except RuntimeProvisionError as exc:
        raise SystemExit(f"fp-tools-runtime: error: {exc}") from exc
    location = activation.distro or str(activation.prefix)
    print(f"{args.component} runtime ready: {location}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
