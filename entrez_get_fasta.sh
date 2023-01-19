#!/usr/bin/env bash
# Take all parameters as a single string
organism=$@
protein="pcrA"
# Replace spaces with +
organism=$(echo "$organism" | sed 's/ /%20/g')
output_file="thermophiles.fasta"

# Store the list of IDs in a variable
ids=$(curl --silent "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=protein&term=$organism+$protein&retmode=json" | jq ".esearchresult.idlist[]" | xargs -I'{}' echo '{}')

# Make the list a comma-separated string
ids=$(echo "$ids" | tr '\n' ',')

# Use jq, iterate over each key and look for the title field, print the value as well as the id
# \n at the end is because we can make it into an array
ids_list=$(curl --silent "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=protein&id=$ids&retmode=json" | jq -r '.result.uids[] as $id | ($id + ";" + .result[$id].title + "\n")')

# List each ID and title, and ask the user to select one, keeping the ID in a variable
# fzf returns id, given "id;titles" string
# selected_id=$(echo "$ids_list" | fzf --preview 'echo {} | cut -d ";" -f 2' | cut -d ";" -f 1)
# Do the same but remove any empty lines
# selected_id=$(echo "$ids_list" | fzf --preview 'echo {} | cut -d ";" -f 2' | cut -d ";" -f 1 | sed '/^$/d')

# Turn ids_list into an array
IFS=$'\n' read -rd '' -a ids_array <<< "$ids_list"
# List each ID and title, and ask the user to select one, keeping the ID in a variable
select id in "${ids_array[@]}"; do
	# Get the ID from the selected line
	id=$(echo "$id" | cut -d ";" -f 1)
	break
done



# Download the sequence in FASTA format 
# check if id is empty
if [ -z "$id" ]; then
	echo "No ID selected"
else
	echo "Downloading sequence $id"
	curl --silent "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&id=$id&rettype=fasta&retmode=text" | tee -a $output_file #output to terminal and file
	echo "Sequence written to $output_file"

	##### Optional
	# Ask for a temperature value 
	read -p "Enter a temperature value: " temp 

	# In the output file, 
	# Get the line number of the last occurence of '>'
	last_seq=$(grep -n ">" $output_file | tail -n 1 | cut -d ":" -f 1)

	# Add the temperature value to the end of the line
	sed -i "${last_seq}s/$/ $temp/" $output_file

fi







# Copy it to clipboard for easy pasting
# cat /tmp/seq | xclip -selection clipboard
