import json
import pandas as pd
from watchdog.observers import Observer
from watchdog.events import FileSystemEventHandler

json_file_path = '/Users/ejm6426/Desktop/rtRBA/rtRBA-main/build_model/input/precursors/import-file-test.json'
excel_file_path = '/Users/ejm6426/Desktop/rtRBA/rtRBA-main/build_model/input/precursors/import-file-test.xlsx'


def convert_json_to_excel():
    with open(json_file_path, 'r') as json_file:
        data = json.load(json_file)

    # Flatten the JSON data and create a DataFrame
    df = pd.json_normalize(data)

    # Write the DataFrame to an Excel file
    df.to_excel(excel_file_path, index=False)


class JSONFileHandler(FileSystemEventHandler):
    def on_modified(self, event):
        if event.src_path == json_file_path:
            print(f"Detected changes in {json_file_path}. Updating Excel file...")
            convert_json_to_excel()


if __name__ == "__main__":
    # Initial conversion from JSON to Excel
    convert_json_to_excel()

    # Set up file system event handler
    event_handler = JSONFileHandler()
    observer = Observer()
    observer.schedule(event_handler, path='path/to/your/json/directory', recursive=False)
    observer.start()

    try:
        print("Watching for changes in the JSON file. Press Ctrl+C to stop.")
        observer.join()
    except KeyboardInterrupt:
        observer.stop()
        observer.join()
