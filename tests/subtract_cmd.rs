use std::io::Write;
use std::process::Command;

use sourmash::signature::{Signature, SigsTrait};

use assert_cmd::prelude::*;
use predicates::str::contains;
use tempfile::{NamedTempFile, TempDir};

#[test]
fn subtract_protein() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = Command::cargo_bin("subtract")?;

    let tmp_dir = TempDir::new()?;

    let mut file = NamedTempFile::new()?;
    writeln!(file, "data/GCA_001593925.1_ASM159392v1_protein.faa.gz.sig")?;
    writeln!(file, "data/GCA_001593935.1_ASM159393v1_protein.faa.gz.sig")?;

    cmd.arg("data/GCA_001593925.1_ASM159392v1_protein.faa.gz.sig")
        .arg(file.path())
        .args(&["-k", "57"])
        .args(&["-s", "100"])
        .args(&["-o", tmp_dir.path().to_str().unwrap()])
        .assert()
        .success();

    let path = tmp_dir
        .path()
        .join("57")
        .join("GCA_001593925.1_ASM159392v1_protein.faa.gz.sig");
    assert!(path.exists());
    let mh = &Signature::from_path(path)?.swap_remove(0).sketches()[0];
    assert_eq!(mh.size(), 0);

    let path = tmp_dir
        .path()
        .join("57")
        .join("GCA_001593935.1_ASM159393v1_protein.faa.gz.sig");
    assert!(path.exists());
    let mh = &Signature::from_path(path)?.swap_remove(0).sketches()[0];
    assert_eq!(mh.size(), 4463);

    Ok(())
}

#[test]
fn subtract_dayhoff() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = Command::cargo_bin("subtract")?;

    let tmp_dir = TempDir::new()?;

    let mut file = NamedTempFile::new()?;
    writeln!(file, "data/bat2-LU__AAACCTGAGCCACGCT.sig")?;

    cmd.arg("data/bat2-LU__AAACCTGAGCCACGCT.sig")
        .arg(file.path())
        .args(&["-k", "42"])
        .args(&["-s", "10"])
        .args(&["-e", "dayhoff"])
        .args(&["-o", tmp_dir.path().to_str().unwrap()])
        .assert()
        .success();

    let path = tmp_dir
        .path()
        .join("42")
        .join("bat2-LU__AAACCTGAGCCACGCT.sig");
    assert!(path.exists());
    let mh = &Signature::from_path(path)?.swap_remove(0).sketches()[0];
    assert_eq!(mh.size(), 0);

    Ok(())
}

#[test]
fn subtract_protein_from_dayhoff() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = Command::cargo_bin("subtract")?;

    let tmp_dir = TempDir::new()?;

    let mut file = NamedTempFile::new()?;
    writeln!(file, "data/GCA_001593935.1_ASM159393v1_protein.faa.gz.sig").unwrap();

    cmd.arg("data/bat2-LU__AAACCTGAGCCACGCT.sig")
        .arg(file.path())
        .args(&["-k", "42"])
        .args(&["-s", "10"])
        .args(&["-e", "dayhoff"])
        .args(&["-o", tmp_dir.path().to_str().unwrap()])
        .assert()
        .failure()
        .stderr(contains(
            r#"Unable to load a sketch from "data/GCA_001593935.1_ASM159393v1_protein.faa.gz.sig"#,
        ));

    Ok(())
}
