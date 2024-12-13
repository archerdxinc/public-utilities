# Brent Lutz, IDT
# This script summarizes key job data from an v6 instance,
# based on a (optionally) defined time period, outputting a CSV.
# Tested and works on 6.2 but also 7.4.2

# REQUIRED: Ensure requests library is installed - pip install requests
# REQUIRED: Ensure env variables API_USERNAME and API_PASSWORD are set first
# REQUIRED: Replace domain url
# OPTIONAL: Specify start_date_str and end_date_str
# Endpoint is specific to v6 "rest_api", which is just "api" in v7
# Modify the page and page_size parameters if desired
# The output CSV will be saved in the current directory where the script is run

import os
import json
import csv
import requests
from requests.auth import HTTPBasicAuth
from datetime import datetime

# Base API domain (replace with the desired domain)
domain = 'https://staging.analysis.archerdx.xyz/'

# Fetch username and password from environment variables
username = os.getenv('API_USERNAME')
password = os.getenv('API_PASSWORD')

# Check if the environment variables are set
if not username or not password:
    raise ValueError("Environment variables API_USERNAME and API_PASSWORD must be set.")

# Base API URLs
base_jobs_url = f'{domain}/rest_api/jobs/'
job_detail_url_template = f'{domain}/rest_api/jobs/{{job_id}}'

# Output file path
output_file = 'filtered_job_details_with_samples.csv'

# User-defined date range (optional)
# Replace with desired date strings in 'YYYY-MM-DD' format, or leave as None
start_date_str = '2024-10-25'  # e.g., '2024-01-01'
end_date_str = '2024-12-10'    # e.g., '2024-12-31'

# Parse dates if provided
start_date = datetime.fromisoformat(start_date_str).date() if start_date_str else None
end_date = datetime.fromisoformat(end_date_str).date() if end_date_str else None

# Initialize variables
job_details = []
page = 1
page_size = 100

# Set up Basic Authentication
auth = HTTPBasicAuth(username, password)

# Fetch paginated jobs
while True:
    # Fetch data from the jobs API for the current page
    api_url = f"{base_jobs_url}?page={page}&page_size={page_size}"
    response = requests.get(api_url, auth=auth)

    if response.status_code != 200:
        print(f"Failed to fetch data from API. Status code: {response.status_code}, Message: {response.text}")
        break

    data = response.json()

    # Extract required fields from the current page's results
    for entry in data['results']:
        job_id = entry['job_id']
        complete_time = entry['job_db_record'].get('complete_time', '')
        
        # Parse complete_time and filter by date range
        complete_date = datetime.fromisoformat(complete_time).date() if complete_time else None
        if start_date and complete_date and complete_date < start_date:
            continue
        if end_date and complete_date and complete_date > end_date:
            continue

        user_email = entry['job_db_record'].get('user_email', '')
        job_status_name = entry['job_db_record'].get('job_status_name', '')
        job_name = entry['job_db_record'].get('job_name', '')
        target_region_name = entry.get('target_region_name', '')

        # Fetch job details to get the samples count
        job_detail_url = job_detail_url_template.format(job_id=job_id)
        job_detail_response = requests.get(job_detail_url, auth=auth)

        if job_detail_response.status_code != 200:
            print(f"Failed to fetch job details for job_id {job_id}. Status code: {job_detail_response.status_code}, Message: {job_detail_response.text}")
            samples_count = None
        else:
            job_detail_data = job_detail_response.json()
            samples_count = len(job_detail_data.get('samples', []))

        job_details.append([
            job_id, complete_date, user_email, job_status_name,
            job_name, target_region_name, samples_count
        ])

    # Check if there is a next page
    if not data.get('next'):
        break

    # Move to the next page
    page += 1

# Write all collected data to a CSV file
with open(output_file, 'w', newline='') as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow(['job_id', 'complete_date', 'user_email', 'job_status_name', 'job_name', 'target_region_name', 'samples_count'])  # Write header
    writer.writerows(job_details)

print(f"Filtered CSV file with job details and sample counts has been saved to {output_file}")
