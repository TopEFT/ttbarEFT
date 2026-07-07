import os
import json

def extract_paths():
    # Get all files in the current directory
    files = os.listdir('.')
    
    # Filter for files ending in _ptSkim.json
    target_files = [f for f in files if f.endswith('_ptSkim.json')]
    
    if not target_files:
        print("No files matching *_ptSkim.json found.")
        return

    print(f"{'FILENAME':<30} | {'PATH VALUE'}")
    print("-" * 80)

    for filename in target_files:
        try:
            with open(filename, 'r') as f:
                data = json.load(f)
                # Use .get() to avoid crashing if "path" is missing in a file
                path_val = data.get('path', 'KEY NOT FOUND')
                print(f"{filename:<30} | {path_val}")
        except Exception as e:
            print(f"Error reading {filename}: {e}")

if __name__ == "__main__":
    extract_paths()
