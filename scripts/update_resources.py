#!/usr/bin/env python3
"""Download, convert, and install ClinVar resources for ANNOVAR."""

from __future__ import annotations

import argparse
import hashlib
import re
import shutil
import subprocess
import sys
import time
from pathlib import Path

import certifi
import requests
import urllib3

import avinput2annovardb


def log(level: str, message: str) -> None:
    stamp = time.strftime("%Y-%m-%d %H:%M:%S")
    print(f"[{level}] {stamp} - {message}")


def fail(message: str) -> None:
    raise SystemExit(f"ERROR: {message}")


def md5(path: Path) -> str:
    digest = hashlib.md5()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def download(url: str, destination: Path) -> None:
    log("INFO", f"Downloading {url}")
    response = requests.get(url, stream=True, timeout=120)
    response.raise_for_status()
    with destination.open("wb") as handle:
        for block in response.iter_content(1024 * 1024):
            if block:
                handle.write(block)


def latest_remote(http: urllib3.PoolManager, url: str) -> tuple[str, str]:
    response = http.request("GET", url)
    if response.status != 200:
        raise RuntimeError(f"ClinVar index returned HTTP {response.status}")
    text = response.data.decode("utf-8")
    releases = re.findall(r'clinvar_(\d+)\.vcf\.gz', text)
    if not releases:
        raise RuntimeError("No ClinVar VCF release found")
    release = max(releases)
    checksum_response = http.request("GET", f"{url}clinvar_{release}.vcf.md5")
    if checksum_response.status != 200:
        raise RuntimeError(f"ClinVar checksum returned HTTP {checksum_response.status}")
    checksum_match = re.search(r"([0-9a-fA-F]{32})", checksum_response.data.decode("utf-8"))
    if not checksum_match:
        raise RuntimeError("ClinVar checksum file has no MD5 value")
    return release, checksum_match.group(1).lower()


def ensure_release(cache: Path, build: str) -> tuple[Path, str]:
    cache.mkdir(parents=True, exist_ok=True)
    url = f"https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_{build}/"
    http = urllib3.PoolManager(cert_reqs="CERT_REQUIRED", ca_certs=certifi.where())
    release, remote_md5 = latest_remote(http, url)
    vcf = cache / f"clinvar_{release}.vcf.gz"
    if not vcf.exists() or md5(vcf) != remote_md5:
        download(f"{url}{vcf.name}", vcf)
        if md5(vcf) != remote_md5:
            vcf.unlink(missing_ok=True)
            raise RuntimeError(f"MD5 verification failed for {vcf}")
    for suffix in (".md5", ".tbi"):
        target = cache / f"{vcf.name}{suffix}"
        if not target.exists():
            download(f"{url}{target.name}", target)
    return vcf, release


def install_clinvar(args: argparse.Namespace) -> Path:
    build_prefix = {"GRCh37": "hg19", "GRCh38": "hg38"}[args.genome_version]
    annovar = Path(args.annovar_path).expanduser().resolve()
    humandb = (Path(args.humandb_path).expanduser().resolve()
               if args.humandb_path else annovar / "humandb")
    humandb.mkdir(parents=True, exist_ok=True)
    convert = annovar / "convert2annovar.pl"
    index = annovar / "index_annovar.pl"
    if not convert.is_file() or not index.is_file():
        fail(f"ANNOVAR scripts not found under {annovar}")

    cache = Path(__file__).resolve().parent / "clinvar" / args.genome_version
    vcf, release = ensure_release(cache, args.genome_version)
    avinput = cache / f"clinvar_{release}.avinput"
    converted = cache / f"clinvar_{release}.txt"
    if not avinput.exists():
        subprocess.run(["perl", str(convert), "-format", "vcf4", "-includeinfo", str(vcf), "-outfile", str(avinput)], check=True)
    if not converted.exists():
        avinput2annovardb.clinvaravinput2annovardb(avinput, avinput2annovardb.DEFAULT_FIELDS, converted)

    name = args.rename or release
    installed = humandb / f"{build_prefix}_clinvar_{name}.txt"
    shutil.move(str(converted), str(installed))
    index_prefix = installed.with_suffix("")
    subprocess.run(["perl", str(index), str(installed), "-outfile", str(index_prefix)], check=True)
    if not any(installed.parent.glob(installed.stem + ".idx")):
        log("WARNING", f"ANNOVAR index output was not found beside {installed}")
    if args.rename:
        (humandb / f"{build_prefix}_clinvar_{args.rename}.ver").write_text(f"ClinVar:{release}\n", encoding="utf-8")
    log("INFO", f"Installed ClinVar database: {installed}")
    return installed


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-d", "--database-type", required=True)
    parser.add_argument("-hp", "--humandb-path")
    parser.add_argument("-g", "--genome-version", choices=("GRCh37", "GRCh38"), required=True)
    parser.add_argument("-a", "--annovar-path", required=True)
    parser.add_argument("-r", "--rename")
    args = parser.parse_args()
    if args.database_type != "clinvar":
        fail("Only the clinvar database is supported")
    try:
        install_clinvar(args)
    except (OSError, RuntimeError, subprocess.CalledProcessError, requests.RequestException) as exc:
        fail(str(exc))


if __name__ == "__main__":
    main()
