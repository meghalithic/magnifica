# Function to convert vector to norm length 
# sum of the squares of the vector is its magnitude
f.normalize_vector <- function(vector) {
    norm_length <- sqrt(sum(vector^2))
    normalized_vector <- vector / norm_length
    return(normalized_vector)
}