use vcpkg;

fn main() {
    let target = std::env::var("CARGO_CFG_TARGET_OS").unwrap();
    if target == "windows" {
        use vcpkg;
        vcpkg::Config::new()
            .find_package("gdal")
            .unwrap();
    }
}
