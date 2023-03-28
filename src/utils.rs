use std::fmt::Display;

pub fn generate_vector_space_delimited<T: Display>(vec: Vec<T>) -> String {
    let mut string = "".to_string();
    for val in vec {
        string.push_str(&format!("{val} "));
    }
    let string = string.trim().to_string();
    string
}
