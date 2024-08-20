import os

def is_unicode(char):
    try:
        char.encode('utf-8')
        return True
    except UnicodeEncodeError:
        return False

def remove_non_unicode(file_path):
    # Read the file's content
    with open(file_path, 'r', encoding='utf-8', errors='ignore') as file:
        content = file.read()

    # Filter out non-Unicode characters
    cleaned_content = ''.join(filter(is_unicode, content))

    # Write the cleaned content back to the file
    with open(file_path, 'w', encoding='utf-8') as file:
        file.write(cleaned_content)

if __name__ == "__main__":
    file_path = "/home/ilianolhin/git/nerpa2/configs/norineText.csv"  # Change to your file's path
    if os.path.exists(file_path):
        remove_non_unicode(file_path)
        print(f"Non-Unicode characters have been removed from {file_path}.")
    else:
        print(f"The file {file_path} does not exist.")
