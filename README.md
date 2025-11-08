# PubMed-MCP

A Model Context Protocol (MCP) server that provides tools for searching PubMed articles using the NCBI Entrez API.

**Author:** Emilio Delgado Muñoz

<a href="https://glama.ai/mcp/servers/@emi-dm/PubMed-MCP">
  <img width="380" height="200" src="https://glama.ai/mcp/servers/@emi-dm/PubMed-MCP/badge" alt="PubMed Server MCP server" />
</a>

## Features

- Search PubMed for articles based on queries
- Retrieve detailed information including title, authors, abstract, journal, and publication date
- Returns results in JSON format
- Configurable maximum number of results

## Architecture

```mermaid
graph TB
  A[User] --> B[MCP Server<br/>pubmed_server.py]
  B --> C[search_pubmed function]
  C --> D[Entrez.esearch<br/>Search in PubMed]
  D --> E[PubMed database<br/>NCBI]
  E --> F[List of PMIDs]
  F --> G[Entrez.efetch<br/>Fetch details]
  G --> E
  G --> H[Article XML records]
  H --> I[Data processing]
  I --> J[Extraction of:<br/>- Title<br/>- Authors<br/>- Abstract<br/>- Journal<br/>- Date]
  J --> K[List of articles<br/>in JSON format]
  K --> L[Response to user]

  subgraph "Dependencies"
    M[BioPython<br/>requirements.txt]
    N[FastMCP<br/>requirements.txt]
  end

  B -.-> M
  B -.-> N

  subgraph "Configuration"
    O[Entrez.email<br/>Configured in code]
  end

  C -.-> O

  style A fill:#e1f5fe
  style L fill:#c8e6c9
  style E fill:#fff3e0
```

## Installation

1. Clone this repository:
   ```bash
   git clone <repository-url>
   cd PubMed-MCP
   ```

2. Install dependencies:
   ```bash
   uv sync
   ```

3. Configure your email in `pubmed_server.py`:
   ```python
   Entrez.email = 'your-email@example.com'  # Replace with your actual email
   ```

## VS Code Configuration

To use this MCP server locally in VS Code, the project includes a pre-configured `.vscode/mcp.json` file. This file tells VS Code how to run the MCP server.

The configuration is already set up to use `uv` for running the server:

```json
{
  "servers": {
    "pubmed-mcp": {
      "command": "uv",
      "args": ["run", "${workspaceFolder}/pubmed_server.py"]
    }
  }
}
```

### Requirements for VS Code Integration

- VS Code with MCP extension support
- `uv` package manager installed
- Python virtual environment set up

### Alternative Configuration

If you prefer to use `pip` instead of `uv`, you can modify the `.vscode/mcp.json` file:

```json
{
  "servers": {
    "pubmed-mcp": {
      "command": "python",
      "args": ["${workspaceFolder}/pubmed_server.py"]
    }
  }
}
```

Make sure your virtual environment is activated when using this configuration.

## Requirements

- Python 3.11+
- BioPython
- FastMCP

## Usage

Run the MCP server:

```bash
python pubmed_server.py
```

The server will start and listen for MCP protocol messages on stdin/stdout.

## Available Tools

### search_pubmed

Searches PubMed for articles matching the given query.

**Parameters:**
- `query` (string): The search query
- `max_results` (integer, optional): Maximum number of results to return (default: 10)
 - `title` (bool, optional): If true (default) search in Title field
 - `abstract` (bool, optional): If true (default) search in Abstract field
 - `keywords` (bool, optional): If true (default) expand search with Author Keywords (`[ot]`) and MeSH Headings (`[mh]`)

**Field logic:**
- `title=True` and `abstract=True` -> query applied as `(your terms)[tiab]`
- Only `title=True` -> `(your terms)[ti]`
- Only `abstract=True` -> `(your terms)[ab]`
- Both false -> no field tag (all fields)
- `keywords=True` -> OR-expanded with `(your terms)[ot] OR (your terms)[mh]`

**Example refined queries:**

```text
query = "breast cancer metastasis"
title=True, abstract=True, keywords=True -> (breast cancer metastasis)[tiab] OR ((breast cancer metastasis)[ot] OR (breast cancer metastasis)[mh])
title=True, abstract=False, keywords=False -> (breast cancer metastasis)[ti]
title=False, abstract=False, keywords=True -> (breast cancer metastasis) OR ((breast cancer metastasis)[ot] OR (breast cancer metastasis)[mh])
```

**Returns:**
A list of article objects containing:

- `pmid`: PubMed ID
- `title`: Article title
- `authors`: List of author names
- `abstract`: Article abstract
- `journal`: Journal name
- `publication_year`: Year of publication
- `publication_month`: Month of publication
- `url`: PubMed URL

## Configuration

Before using the tool, you must set your email address in the `Entrez.email` variable. This is required by NCBI's Entrez API.

## License

This project is open source. Please check the license file for details.

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.