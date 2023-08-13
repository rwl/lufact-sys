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
    println!("cargo:rustc-link-lib=gfortran");

    let output = Command::new("make")
        .args(&["-C", vendor_path.to_str().unwrap(), "all"])
        .output()
        .expect("could not spawn `make`");
    assert!(output.status.success(), "{:?}", output);
}
