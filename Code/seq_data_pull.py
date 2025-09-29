# from Bio import Entrez
# import os
# import sys
# import time

# # Always tell NCBI who you are
# Entrez.email = "your.email@example.com"

# def get_sra_ids_from_samn(samn_ids):
#     """
#     Fetches SRA run IDs (SRR) for a list of BioSample accession numbers (SAMN).
#     """
#     sra_ids = []
    
#     # Batch search for efficiency
#     id_string = " OR ".join(samn_ids)
    
#     try:
#         # Search the BioSample database first
#         handle = Entrez.esearch(db="biosample", term=id_string, retmax=len(samn_ids))
#         biosample_record = Entrez.read(handle)
#         handle.close()
        
#         biosample_ids = biosample_record["IdList"]
#         print(f"Found {len(biosample_ids)} BioSample IDs.")

#         if not biosample_ids:
#             return []

#         # Now, link the BioSample IDs to the SRA database to get SRR IDs
#         handle = Entrez.elink(dbfrom="biosample", db="sra", id=",".join(biosample_ids))
#         link_record = Entrez.read(handle)
#         handle.close()

#         for link_set in link_record:
#             if "LinkSetDb" in link_set and link_set["LinkSetDb"]:
#                 for link_info in link_set["LinkSetDb"]:
#                     if "Link" in link_info:
#                         for link in link_info["Link"]:
#                             sra_ids.append(link["Id"])

#     except Exception as e:
#         print(f"An error occurred: {e}", file=sys.stderr)
    
#     return sra_ids

# def write_commands_to_file(sra_ids, filename):
#     """
#     Writes prefetch and fastq-dump commands to a shell script file.
#     """
#     # Use "w" mode to create or overwrite the file
#     with open(filename, "w") as f:
#         f.write("#!/bin/bash\n")
#         f.write("# This script downloads and converts SRA files to FASTQ.\n")
#         f.write("# Please make sure the SRA Toolkit is installed and in your PATH.\n\n")
        
#         # Write prefetch command
#         f.write(f"prefetch --max-size u {' '.join(sra_ids)}\n")
        
#         # Write fastq-dump command
#         f.write(f"fastq-dump --split-files --gzip {' '.join(sra_ids)}\n")
    
#     print(f"Commands written to {filename}. You can run it from your terminal using 'bash {filename}'")

# def main():
#     # Generate lists of accession numbers
#     gyrB_samn_accessions = [f"SAMN{i}" for i in range(10964863, 10965439)]
#     s16S_samn_accessions = [f"SAMN{i}" for i in range(10970131, 10970567)]

#     # Fetch SRA IDs for both amplicon types
#     print("Fetching SRA accessions for gyrB data...")
#     gyrB_sra_ids = get_sra_ids_from_samn(gyrB_samn_accessions)
#     print(f"Found {len(gyrB_sra_ids)} gyrB SRA accessions.")
    
#     print("Fetching SRA accessions for 16S data...")
#     s16S_sra_ids = get_sra_ids_from_samn(s16S_samn_accessions)
#     print(f"Found {len(s16S_sra_ids)} 16S SRA accessions.")
    
#     # Write the commands to a file
#     write_commands_to_file(gyrB_sra_ids, "download_gyrB_data.sh")
#     write_commands_to_file(s16S_sra_ids, "download_16S_data.sh")

# if __name__ == "__main__":
#     main()

# from Bio import Entrez
# import os
# import sys
# import time

# # Always tell NCBI who you are
# Entrez.email = "your.email@example.com"

# def get_sra_ids_from_samn(samn_ids):
#     """
#     Fetches SRA run IDs (SRR) for a list of BioSample accession numbers (SAMN).
#     """
#     sra_ids = []
    
#     try:
#         # Step 1: Search the SRA database for the BioSample accessions.
#         # This is more direct than searching BioSample first.
#         term_query = " OR ".join(samn_ids)
#         handle = Entrez.esearch(db="sra", term=term_query, retmax=100000)
#         sra_record = Entrez.read(handle)
#         handle.close()
        
#         sra_ids = sra_record["IdList"]
#         print(f"Found {len(sra_ids)} SRA IDs.")

#     except Exception as e:
#         print(f"An error occurred: {e}", file=sys.stderr)
    
#     return sra_ids

# def write_commands_to_file(sra_ids, filename):
#     """
#     Writes prefetch and fastq-dump commands to a shell script file.
#     """
#     if not sra_ids:
#         print(f"No SRA IDs to write for {filename}.")
#         return

#     # Use "w" mode to create or overwrite the file
#     with open(filename, "w") as f:
#         f.write("#!/bin/bash\n")
#         f.write("# This script downloads and converts SRA files to FASTQ.\n")
#         f.write("# Please make sure the SRA Toolkit is installed and in your PATH.\n\n")
        
#         # Write prefetch command
#         f.write(f"prefetch --max-size u {' '.join(sra_ids)}\n")
        
#         # Write fastq-dump command
#         f.write(f"fastq-dump --split-files --gzip {' '.join(sra_ids)}\n")
    
#     print(f"Commands written to {filename}. You can run it from your terminal using 'bash {filename}'")

# def main():
#     # Generate lists of accession numbers
#     # Note: These ranges may not contain actual metagenomic data.
#     # The code's logic is what has been fixed.
#     # gyrB_samn_accessions = [f"SAMN{i}" for i in range(10964863, 10965439)] #filter this out
#     s16S_samn_accessions = [f"SAMN{i}" for i in range(10970131, 10970567)]

#     # Fetch SRA IDs for both amplicon types
#     # print("Fetching SRA accessions for gyrB data...")
#     # gyrB_sra_ids = get_sra_ids_from_samn(gyrB_samn_accessions)
    
#     print("Fetching SRA accessions for 16S data...")
#     s16S_sra_ids = get_sra_ids_from_samn(s16S_samn_accessions)
    
#     # Write the commands to a file
#     # write_commands_to_file(gyrB_sra_ids, "download_gyrB_data.sh")
#     write_commands_to_file(s16S_sra_ids, "download_16S_data.sh")

# if __name__ == "__main__":
#     main()




##############################################################################

# from Bio import Entrez
# import os
# import sys
# import time

# # Always tell NCBI who you are
# # Replace with your actual email
# Entrez.email = "jedshakarji@gmail.com"

# def get_sra_ids_from_samn(samn_ids, batch_size=1000):
#     """
#     Fetches SRA run IDs (SRR) for a list of BioSample accession numbers (SAMN)
#     in batches to handle large lists efficiently.
#     """
#     sra_ids = []
#     total_ids = len(samn_ids)
#     print(f"Starting search for {total_ids} BioSample accessions.")
    
#     # Process in batches to avoid overwhelming the Entrez server
#     for i in range(0, total_ids, batch_size):
#         batch_samn_ids = samn_ids[i:i + batch_size]
#         id_string = " OR ".join(batch_samn_ids)
        
#         try:
#             # Search the SRA database directly for the BioSample accessions.
#             # Entrez.esearch can handle different accession types in the term query.
#             handle = Entrez.esearch(db="sra", term=id_string, retmax=batch_size, idtype="acc")
#             sra_record = Entrez.read(handle)
#             handle.close()
            
#             # The 'IdList' from esearch contains the SRA accession numbers
#             sra_ids.extend(sra_record["IdList"])
#             print(f"  Processed batch {i // batch_size + 1}: Found {len(sra_record['IdList'])} SRA IDs.")
            
#             # Be a good neighbor and don't hit the server too hard
#             time.sleep(0.5)

#         except Exception as e:
#             print(f"An error occurred in a batch: {e}", file=sys.stderr)
            
#     return sra_ids

# def write_commands_to_file(sra_ids, filename):
#     """
#     Writes prefetch and fastq-dump commands to a shell script file.
#     """
#     if not sra_ids:
#         print(f"No SRA IDs to write for {filename}.")
#         return

#     # Use "w" mode to create or overwrite the file
#     with open(filename, "w") as f:
#         f.write("#!/bin/bash\n")
#         f.write("# This script downloads and converts SRA files to FASTQ.\n")
#         f.write("# Please make sure the SRA Toolkit is installed and in your PATH.\n\n")
        
#         # Write prefetch command with a reasonable batch size
#         # This prevents the command line from getting too long
#         batch_size = 500
#         for i in range(0, len(sra_ids), batch_size):
#             batch_sra_ids = sra_ids[i:i + batch_size]
#             f.write(f"echo \"Prefetching batch {i // batch_size + 1}...\"\n")
#             f.write(f"prefetch --max-size u {' '.join(batch_sra_ids)}\n\n")

#         # Write fastq-dump command
#         f.write("echo \"Converting .sra files to .fastq.gz...\"\n")
#         f.write(f"fastq-dump --split-files --gzip --skip-technical --clip {' '.join(sra_ids)}\n")
    
#     print(f"Commands written to {filename}. You can run it from your terminal using 'bash {filename}'")

# def main():
#     # Define the ranges of BioSample accession numbers
#     # These are illustrative ranges from the original user prompt.
#     # gyrB_samn_accessions = [f"SAMN{i}" for i in range(10964863, 10965439)]
#     s16S_samn_accessions = [f"SAMN{i}" for i in range(10970131, 10970567)]

#     # Fetch SRA IDs for both amplicon types
#     print("Fetching SRA accessions for gyrB data...")
#     # gyrB_sra_ids = get_sra_ids_from_samn(gyrB_samn_accessions)
#     # print(f"Found a total of {len(gyrB_sra_ids)} gyrB SRA accessions.")
    
#     print("\nFetching SRA accessions for 16S data...")
#     s16S_sra_ids = get_sra_ids_from_samn(s16S_samn_accessions)
#     print(f"Found a total of {len(s16S_sra_ids)} 16S SRA accessions.")
    
#     # Write the commands to a file
#     # write_commands_to_file(gyrB_sra_ids, "download_gyrB_data.sh")
#     write_commands_to_file(s16S_sra_ids, "download_16S_data.sh")

# if __name__ == "__main__":
#     main()



from Bio import Entrez
import os
import sys
import time
import csv

# Always tell NCBI who you are
# Replace with your actual email
Entrez.email = "your.email@example.com"

def get_sra_ids_from_samn(samn_ids, batch_size=200):
    """
    Finds BioSample IDs, links to SRA, and fetches SRR accessions using the
    robust efetch 'runinfo' method.
    """
    srr_accessions = []
    total_ids = len(samn_ids)
    print(f"Starting search for {total_ids} BioSample accessions.")
    
    for i in range(0, total_ids, batch_size):
        batch_samn_ids = samn_ids[i:i + batch_size]
        id_string = " OR ".join(batch_samn_ids)
        
        try:
            # --- Step 1: Search the BioSample database ---
            handle = Entrez.esearch(db="biosample", term=id_string, retmax=batch_size)
            biosample_record = Entrez.read(handle)
            handle.close()
            biosample_uids = biosample_record["IdList"]

            if not biosample_uids:
                print(f"  Processed batch {i // batch_size + 1}: No valid BioSamples found.")
                continue

            # --- Step 2: Link from BioSample to the SRA database ---
            handle = Entrez.elink(dbfrom="biosample", db="sra", id=",".join(biosample_uids))
            link_records = Entrez.read(handle)
            handle.close()
            
            sra_uids = []
            for record in link_records:
                if 'LinkSetDb' in record and record['LinkSetDb']:
                    for link_db in record['LinkSetDb']:
                        if link_db.get('DbTo') == 'sra':
                            for link in link_db['Link']:
                                sra_uids.append(link['Id'])

            if not sra_uids:
                print(f"  Processed batch {i // batch_size + 1}: Found BioSamples but no linked SRA records.")
                continue

            # --- Step 3: Use efetch to get the RunInfo table (CSV) ---
            handle = Entrez.efetch(db="sra", id=",".join(sra_uids), rettype='runinfo', retmode='text')
            
            # *** THIS IS THE FIX: Decode bytes into a string before parsing ***
            csv_data = handle.read().decode('utf-8').strip().splitlines()
            runinfo_reader = csv.reader(csv_data)
            
            header = next(runinfo_reader)
            try:
                run_col_index = header.index('Run')
            except ValueError:
                print(f"  Processed batch {i // batch_size + 1}: Could not find 'Run' column in RunInfo.")
                continue
                
            batch_srr_ids = [row[run_col_index] for row in runinfo_reader]
            srr_accessions.extend(batch_srr_ids)
            print(f"  Processed batch {i // batch_size + 1}: Found {len(batch_srr_ids)} SRR accessions.")
            
            time.sleep(0.5)

        except Exception as e:
            print(f"An error occurred in a batch: {e}", file=sys.stderr)
            
    return srr_accessions

def write_commands_to_file(sra_ids, filename):
    """
    Writes prefetch and fastq-dump commands to a shell script file.
    """
    if not sra_ids:
        print(f"No SRA IDs to write for {filename}.")
        return

    with open(filename, "w") as f:
        f.write("#!/bin/bash\n")
        f.write("# This script downloads and converts SRA files to FASTQ.\n\n")
        
        batch_size = 500
        for i in range(0, len(sra_ids), batch_size):
            batch_sra_ids = sra_ids[i:i + batch_size]
            f.write(f"echo \"Prefetching batch {i // batch_size + 1}...\"\n")
            f.write(f"prefetch --max-size u {' '.join(batch_sra_ids)}\n\n")

        f.write("echo \"\nConverting .sra files to .fastq.gz...\"\n")
        f.write(f"fastq-dump --split-files --gzip --skip-technical --clip {' '.join(sra_ids)}\n")
    
    print(f"Commands written to {filename}. You can run it from your terminal using 'bash {filename}'")

def main():
    gyrB_samn_accessions = [f"SAMN{i}" for i in range(10964863, 10965439)]
    s16S_samn_accessions = [f"SAMN{i}" for i in range(10970131, 10970567)]

    print("--- Fetching SRA accessions for gyrB data ---")
    gyrB_sra_ids = get_sra_ids_from_samn(gyrB_samn_accessions)
    print(f"Found a total of {len(gyrB_sra_ids)} gyrB SRA accessions.\n")
    
    print("--- Fetching SRA accessions for 16S data ---")
    s16S_sra_ids = get_sra_ids_from_samn(s16S_samn_accessions)
    print(f"Found a total of {len(s16S_sra_ids)} 16S SRA accessions.\n")
    
    write_commands_to_file(gyrB_sra_ids, "download_gyrB_data.sh")
    write_commands_to_file(s16S_sra_ids, "download_16S_data.sh")

if __name__ == "__main__":
    main()