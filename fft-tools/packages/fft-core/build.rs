fn main() {
    println!("cargo:rustc-link-lib=dylib=blas");
    println!("cargo:rustc-link-lib=dylib=lapack");
}
