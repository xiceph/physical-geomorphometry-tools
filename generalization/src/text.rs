pub fn bold(text: String) -> String {
    let modified_text = text.replace("\x1b[0m", "\x1b[0m\x1b[1m");
    format!("\x1b[1m{}\x1b[0m", modified_text)
}

pub fn error(text: String) -> String {
    let modified_text = text.replace("\x1b[0m", "\x1b[0m\x1b[31m");
    format!("\x1b[31m{}\x1b[0m", modified_text)
}

pub fn warning(text: String) -> String {
    let modified_text = text.replace("\x1b[0m", "\x1b[0m\x1b[35;93m");
    format!("\x1b[35;93m{}\x1b[0m", modified_text)
}

pub fn success(text: String) -> String {
    let modified_text = text.replace("\x1b[0m", "\x1b[0m\x1b[92m");
    format!("\x1b[92m{}\x1b[0m", modified_text)
}

pub fn highlight(text: String) -> String {
    let modified_text = text.replace("\x1b[0m", "\x1b[0m\x1b[94m");
    format!("\x1b[94m{}\x1b[0m", modified_text)
}
