SumStats_df <- function(data, vars = names(data), binary = c()) {
  rows <- vars
  cols <- c("Mean (SD, range)", "n")
  out <- matrix("", nrow = length(rows), ncol  = length(cols),
                dimnames = list(rows, cols))
  
  ifelse(length(binary) > 0, rows1 <- rows[-binary], rows1 <- rows)
  for(r in which(rows %in% rows1)) {
    vec <- data %>% pull(rows[r])
    out[rows[r], "Mean (SD, range)"] <- vec %>% mean(na.rm = T) %>% round(digits = 2) %>%
      str_c(" (", vec %>% sd(na.rm = T) %>% round(digits = 2), ", ",
            vec %>% min(na.rm = T) %>% round(digits = 2), "-",
            vec %>% max(na.rm = T) %>% round(digits = 2), ")")
    out[rows[r], "n"] <- sum(!is.na(vec))
  }
  
  if(length(binary) > 0) {
    rows1 <- rows[binary]
    for(r in which(rows %in% rows1)) {
      vec <- data %>% pull(rows[r])
      out[rows[r], "Mean (SD, range)"] <- vec %>%
        (function(x) round(sum(x, na.rm = T) / sum(!is.na(x)) * 100)) %>%
        str_c("%")
      out[rows[r], "n"] <- sum(!is.na(vec))
    }
  }
  
  return(out)
}
