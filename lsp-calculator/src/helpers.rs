use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

// Loads a configuration value from `app.config` for a given key.
pub fn load_config_value(key_to_find: &str) -> Option<String> {
    let config_path = Path::new("app.config");
    if !config_path.exists() {
        return None;
    }

    let file = File::open(config_path).ok()?;
    let reader = BufReader::new(file);

    for line in reader.lines() {
        let line = line.ok()?;
        if let Some((key, value)) = line.split_once('=') {
            if key.trim() == key_to_find {
                return Some(value.trim().to_string());
            }
        }
    }
    None
}
