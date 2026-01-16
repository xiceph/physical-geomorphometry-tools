use console::{style, Emoji};

pub static CHECK: Emoji<'static, 'static> = Emoji("✓", "+");
pub static CROSS: Emoji<'static, 'static> = Emoji("✗", "x");

pub fn check_icon() -> String {
    style(format!("{}", CHECK)).green().to_string()
}

/// Returns a styled string for an arrow icon, with a fallback for non-interactive terminals.
pub fn arrow_icon() -> String {
    if console::user_attended() {
        style("→").cyan().to_string()
    } else {
        style("->").cyan().to_string()
    }
}

pub fn bold<T: AsRef<str>>(text: T) -> String {
    style(text.as_ref()).bold().to_string()
}

pub fn error_icon() -> String {
    style(format!("{}", CROSS)).red().to_string()
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
