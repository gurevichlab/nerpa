import requests
from pathlib import Path
import zipfile
from src.pipeline.logger import NerpaLogger
from bs4 import BeautifulSoup


def get_zip_filename(job_id: str) -> str:
    # Get the .zip filename for the given job ID by scraping the upload directory page.
    # WARNING: this was generated with ChatGPT, I didn't check what it does.
    url = f"https://antismash.secondarymetabolites.org/upload/{job_id}/"

    # Send a request to get the HTML content of the page
    response = requests.get(url)
    response.raise_for_status()

    # Parse the page content
    soup = BeautifulSoup(response.text, 'html.parser')

    # Find all anchor tags with href attributes
    for link in soup.find_all('a', href=True):
        href = link['href']
        if href.endswith(".zip"):
            return href  # Return the first .zip file found

    raise ValueError(f"No .zip file found for job {job_id}")

    
def download_results_zip(job_id: str,
                         output_dir: Path,
                         log: NerpaLogger) -> Path:
    """
    Download the antiSMASH result files for a given job ID once it's complete.

    Parameters:
    job_id (str): The job ID for which results are to be downloaded.
    output_dir (Path): The directory where the results will be saved.
    returns the Path to the downloaded archive
    """
    zip_file_name = get_zip_filename(job_id)
    download_url = f"https://antismash.secondarymetabolites.org/upload/{job_id}/{zip_file_name}"
    log.info(f'Downloading antiSMASH results from {download_url}')

    # Send a request to download the result as a zip file
    response = requests.get(download_url)

    if response.status_code == 200:
        downloaded_zip_filename = output_dir / f"{job_id}.zip"

        # Save the content as a zip file
        with downloaded_zip_filename.open("wb") as file:
            file.write(response.content)

        log.info(f"Results downloaded successfully to {downloaded_zip_filename}")
    else:
        raise ValueError(f"Error: Unable to download results. HTTP Status code: {response.status_code}")

    return downloaded_zip_filename


def download_antismash_results(job_id: str,
                               output_dir: Path,
                               log: NerpaLogger) -> Path:
    """
    Main function to manage the job status checking and downloading of antiSMASH results.

    Parameters:
    job_id (str): The job ID to check the status and download the results for.
    output_dir (Path): The directory where the results will be saved.
    returns path to the folder with the downloaded results
    """

    results_zip = download_results_zip(job_id, output_dir, log)
    results_dir = output_dir / job_id
    results_dir.mkdir()
    with zipfile.ZipFile(results_zip, 'r') as zip_ref:
        # Extract all contents to the specified directory
        zip_ref.extractall(results_dir)
    results_zip.unlink()

    return output_dir


def is_antismash_results_dir(p: Path) -> bool:
    try:
        html_content = BeautifulSoup((p / 'index.html').read_text(), 'html.parser')
        return 'antiSMASH' in html_content.title.string
    except:
        return False


if __name__ == "__main__":
    # Example Job ID (replace with actual job ID)
    job_id = "bacteria-eaecebfb-0eb6-4db1-ac25-5ef2451b4b4a"
    # Run the main function
    download_antismash_results(job_id,  Path('.'), NerpaLogger())