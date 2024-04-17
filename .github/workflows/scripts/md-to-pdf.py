# import httpx
import gh_md_to_html
# from xhtml2pdf import pisa

# Combine markdown files to a single string
# with open("README.md") as mdfile:
#     body = mdfile.read()

# Create table of contents

# Convert to HTML using GitHub API
# response = httpx.post(
#     "https://api.github.com/markdown",
#     json={
#         "mode": "gfm",
#         "text": body
#     }
# )
# print("GitHub API response code:",response.status_code)
# if response.status_code != 200:
#     print("Exit due to wrong API response.")
#     exit(1)

gh_md_to_html.main("-t file -c gh-md-css ../USER_MANUAL.md")

# Add surrounding html for pdf conversion
# (also add that pdf was created using httpx and xhtml2pdf Python modules)


# with open("VHoppAS-documentation.pdf",'w') as pdffile:
#     pdffile.write(htmldoc)