import os
from operator import itemgetter

def convert_hog2og_n0(new_og_dict, header_str, output_path):
    with open(output_path, "w") as writer:
        writer.write(header_str + "\n")
        for new_og_name, ogs in new_og_dict.items():
            writer.write(output_path,"\t".join([new_og_name] + ogs))

def update_hogs(input_path):
    sorted_matrix, header_str = read_old_hogs(input_path) 
    name_dictionary = {}      
    new_og_dict = {}
    # For each line in sorted HOG order replace HOG name with index OG name (based on length of enumerate so HOG.N0 + 0*x + number)
    for pos, line in enumerate(sorted_matrix):
        new_og_name = "OG%07d" % pos
        key = line[2]
        if key in name_dictionary:
            additional = [new_og_name] + [line[1]] + [line[3]]           
            new_value = name_dictionary[key]
            new_value.append(additional)
            name_dictionary.update({key:new_value})
        else:
            name_dictionary[key] = [[new_og_name] + [line[1]] + [line[3]]]
        new_og_dict[new_og_name] = line[4:]

    return name_dictionary, new_og_dict, header_str

def read_old_hogs(input_path):
    #holds lines to write to new output file
    matrix = []
    with open(input_path) as input_file:
        for line in input_file:
            # for each line split into cells (Og, species 1, species 2 etc).
            line_split = line.split("\t")
            # Count number of genes in each line as a count.
            count = len(list(filter(None, (",".join(line_split[3:]).split(",")))))
            # Add line with gene count to matrix variable.
            matrix.append([count] + line_split)

    header = matrix[0]
    matrix.pop(0)
    header.remove("Gene Tree Parent Clade")
    header_str = "\t".join(["Orthogroups"] + header[3:])
        # sort HOG lines by number of genes key
    sorted_matrix =sorted(matrix, key=itemgetter(0), reverse=True)

    return sorted_matrix, header_str

def update_all_hog_files(input_folder, output_folder, name_dictionary):
        # Get all files in the input folder, filtering out N0.tsv
    for file_name in os.listdir(input_folder):
        if file_name.endswith(".tsv") and file_name != "N0.tsv":
            input_file = os.path.join(input_folder, file_name)
            output_file = os.path.join(output_folder, file_name)
            # Call the function for each file
            convert_all_hog_files(input_file, output_file, name_dictionary)  
    

# This converts the N#.tsv so that each N#.HOG references the updated OG id from the new OG file created in convertHOGtoOG
# Each line is check to convert the OG that the HOG is made from into the correct ID.
def convert_all_hog_files(input_path, output_path, name_dictionary):

    with open(input_path) as reader, open(output_path, "w") as writer:
        for line in reader:
            line_split = line.split("\t")
            # excludes header. 
            if line_split[0] != "HOG":
                    ### for each line get the old OG name and find the new set of OG names (N0.HOG)                
                    old_OG_name = line_split[1]
                    key_value_pair = name_dictionary[old_OG_name]

                    ## exception for rows containing 2 genes with no node data...
                    if line_split[2] == "-":
                        node = 0
                    ## Get node id for that HOG this allows us to tell which new OG id belongs to which HOG
                    else:
                        node = int(line_split[2][1:])

                    ## Sort node id in accessending order..
                    sorted_key_value_pair = sorted(key_value_pair, key=itemgetter(1))
                    internal = 0
                    
                    ## Iterate through nodes until the new HOG node is greater than one and less than another... this is the new OG id - i.e. the HOG that this N.HOG comes from
                    ## Set internal as last value incase condition not met until last..
                    internal = len(sorted_key_value_pair) - 1
                    for position,pair in enumerate(sorted_key_value_pair):
                        if pair[-1] == "-":
                            pair_node=0
                        else:
                            pair_node = int(pair[-1][1:])

                        if node < pair_node:
                            internal = position - 1
                            break
                    if internal == -1:
                        internal = 0
                            
                    ## Update that line to replace old OG id with that of the new OG ID from the N0 file.
                    line_split[1] = sorted_key_value_pair[internal][0]
                    writer.write(output_path,"\t".join(line_split))
            # ensures header is written.
            else:
                writer.write(output_path,"\t".join(line_split))  

