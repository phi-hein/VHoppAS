import httpx

with open("README.md") as mdfile:
    body = mdfile.read()

response = httpx.post(
    "https://api.github.com/markdown",
    json={
        "mode": "gfm",
        "text": body
    }
)

print("GitHub API response code:",response.status_code)
if response.status_code != 200:
    print("Exit due to wrong API response.")
    exit(1)

with open(".github/workflows/scripts/VHoppAS-documentation.pdf",'w') as pdffile:
    pdffile.write(response.text)