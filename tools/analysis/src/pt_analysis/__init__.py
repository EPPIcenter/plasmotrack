"""Reader for plasmotrack Run stores and Analysis manifests."""

from importlib.metadata import PackageNotFoundError, version

try:
    __version__ = version("pt-analysis")
except PackageNotFoundError:  # imported from the source tree without being installed
    __version__ = "unknown"
