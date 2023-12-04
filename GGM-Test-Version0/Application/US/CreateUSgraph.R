# Mapping between two-letter abbreviations and full state names
state_mapping <- c(
  "AL" = "Alabama", "AK" = "Alaska", "AZ" = "Arizona", "AR" = "Arkansas", "CA" = "California",
  "CO" = "Colorado", "CT" = "Connecticut", "DE" = "Delaware", "DC" = "District of Columbia",
  "FL" = "Florida", "GA" = "Georgia", "HI" = "Hawaii", "ID" = "Idaho", "IL" = "Illinois",
  "IN" = "Indiana", "IA" = "Iowa", "KS" = "Kansas", "KY" = "Kentucky", "LA" = "Louisiana",
  "ME" = "Maine", "MD" = "Maryland", "MA" = "Massachusetts", "MI" = "Michigan", "MN" = "Minnesota",
  "MS" = "Mississippi", "MO" = "Missouri", "MT" = "Montana", "NE" = "Nebraska", "NV" = "Nevada",
  "NH" = "New Hampshire", "NJ" = "New Jersey", "NM" = "New Mexico", "NY" = "New York",
  "NC" = "North Carolina", "ND" = "North Dakota", "OH" = "Ohio", "OK" = "Oklahoma", "OR" = "Oregon",
  "PA" = "Pennsylvania", "RI" = "Rhode Island", "SC" = "South Carolina", "SD" = "South Dakota",
  "TN" = "Tennessee", "TX" = "Texas", "UT" = "Utah", "VT" = "Vermont", "VA" = "Virginia",
  "WA" = "Washington", "WV" = "West Virginia", "WI" = "Wisconsin", "WY" = "Wyoming"
)

# Read the file
file_path <- '/Users/stahd/Downloads/state_neighbors.txt'
lines <- readLines(file_path)

# Initialize an empty adjacency matrix
n_states <- length(state_mapping)
adj_matrix <- matrix(0, nrow = n_states, ncol = n_states, 
                     dimnames = list((state_mapping), (state_mapping)))

# Populate the adjacency matrix
for(line in lines) {
  tokens <- unlist(strsplit(line, " "))
  state <- state_mapping[tokens[1]]
  neighbors <- state_mapping[tokens[-1]]
  
  for(neighbor in neighbors) {
    adj_matrix[state, neighbor] <- 1
    adj_matrix[neighbor, state] <- 1  # Assuming undirected graph
  }
}

# Display the adjacency matrix
print(adj_matrix)
graph0=adj_matrix
saveRDS(graph0,file = 'us_stat_adj.rds')
