"""Verified network access shared by managed fp-tools downloads."""

from __future__ import annotations

import os
import ssl
import urllib.parse
import urllib.request
from pathlib import Path
from typing import Any

import certifi


CA_ENVIRONMENT_VARIABLES = (
    "FP_TOOLS_CA_BUNDLE",
    "SSL_CERT_FILE",
    "REQUESTS_CA_BUNDLE",
    "CURL_CA_BUNDLE",
)


def verified_ca_bundle() -> Path:
    """Return a valid user-selected CA bundle or the packaged certifi bundle."""

    for variable in CA_ENVIRONMENT_VARIABLES:
        configured = os.environ.get(variable)
        if configured:
            path = Path(configured).expanduser()
            if path.is_file():
                return path.resolve()
    return Path(certifi.where()).resolve()


def verified_ssl_context() -> ssl.SSLContext:
    """Build an SSL context that always requires server-certificate validation."""

    context = ssl.create_default_context(cafile=str(verified_ca_bundle()))
    context.check_hostname = True
    context.verify_mode = ssl.CERT_REQUIRED
    return context


def verified_urlopen(
    request: str | urllib.request.Request,
    *,
    timeout: float,
) -> Any:
    """Open a URL with the explicit fp-tools verification context."""

    return urllib.request.urlopen(
        request,
        timeout=timeout,
        context=verified_ssl_context(),
    )


def is_certificate_failure(error: BaseException) -> bool:
    """Return whether an exception chain represents certificate verification."""

    current: BaseException | None = error
    visited: set[int] = set()
    while current is not None and id(current) not in visited:
        visited.add(id(current))
        if isinstance(current, ssl.SSLCertVerificationError):
            return True
        if "CERTIFICATE_VERIFY_FAILED" in str(current).upper():
            return True
        reason = getattr(current, "reason", None)
        current = reason if isinstance(reason, BaseException) else current.__cause__
    return False


def network_error_message(
    url: str,
    error: BaseException,
    *,
    action: str,
) -> str:
    """Describe a network failure without echoing URL queries or credentials."""

    host = urllib.parse.urlsplit(url).hostname or "the remote server"
    if is_certificate_failure(error):
        return (
            f"TLS certificate verification failed while {action} from {host}. "
            "Install a trusted CA bundle or set FP_TOOLS_CA_BUNDLE or "
            "SSL_CERT_FILE to a valid PEM file."
        )
    reason = getattr(error, "reason", None)
    detail = str(reason or error).strip() or error.__class__.__name__
    return f"Network error while {action} from {host}: {detail}"
