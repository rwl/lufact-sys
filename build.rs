use std::env;
use std::path::PathBuf;
use std::process::Command;

fn main() {
    let out_dir = env::var("OUT_DIR").unwrap();

    let vendor_path = PathBuf::from("vendor")
        .canonicalize()
        .expect("must be able to canonicalize `vendor` path");

    println!("cargo:rustc-link-search={}", out_dir);

    println!("cargo:rustc-link-lib=gp");

    let out = Command::new("make")
        .args(&["-C", vendor_path.to_str().unwrap(), "clean", "install"])
        .output()
        .expect("could not spawn `make`");
    assert!(out.status.success(), "{:?}", out);

    let out = Command::new("make")
        .args(&["-C", vendor_path.to_str().unwrap(), "purge"])
        .output()
        .expect("could not spawn `make`");
    assert!(out.status.success(), "{:?}", out);
}
