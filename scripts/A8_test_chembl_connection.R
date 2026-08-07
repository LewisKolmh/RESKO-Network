library(jsonlite)

url <- "https://www.ebi.ac.uk/chembl/api/data/target/search.json?q=EEF1A1"

result <- fromJSON(url)

str(result)