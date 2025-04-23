use console::style;

pub fn bold<T: AsRef<str>>(text: T) -> String {
    style(text.as_ref()).bold().to_string()
}

pub fn underline<T: AsRef<str>>(text: T) -> String {
    style(text.as_ref()).underlined().to_string()
}

pub fn error<T: AsRef<str>>(text: T) -> String {
    style(text.as_ref()).red().to_string()
}

pub fn warning<T: AsRef<str>>(text: T) -> String {
    style(text.as_ref()).color256(214).bold().to_string()
}

pub fn success<T: AsRef<str>>(text: T) -> String {
    style(text.as_ref()).green().to_string()
}

pub fn highlight<T: AsRef<str>>(text: T) -> String {
    style(text.as_ref()).blue().bold().to_string()
}

pub fn light<T: AsRef<str>>(text: T) -> String {
    style(text.as_ref()).color256(245).to_string()
}

#[cfg(windows)]
pub const ARROW: &str = ">";

#[cfg(not(windows))]
pub const ARROW: &str = "▶";

#[cfg(windows)]
pub const CHECK: &str = "v";

#[cfg(not(windows))]
pub const CHECK: &str = "✓";
