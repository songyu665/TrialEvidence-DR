# Public release checklist

## GitHub

1. Create a public repository named `TrialEvidence-DR`.
2. Upload the contents of this directory, including `LICENSE`, `CITATION.cff`, and `README.md`.
3. Do not upload credentials, raw-data archives, internal logs, or private session files.
4. Verify that the repository is public and that the default branch contains the intended release.
5. Add the repository URL to the manuscript and release metadata after it exists.

## Zenodo

1. Link the GitHub repository to Zenodo or upload a fixed release archive.
2. Review `zenodo.json` and replace the collective author placeholder with the final author list if required.
3. Create a versioned release and archive it through Zenodo.
4. Obtain the real DOI from Zenodo and add it to `CITATION.cff`, the README, and the manuscript.
5. Test the public download and reproduction instructions in a clean environment.

## Current status

- MIT license: published
- GitHub repository: https://github.com/songyu665/TrialEvidence-DR
- Version tag and GitHub Release: `v0.1.0`
- Zenodo archive and DOI: pending
- Public source upload: completed
