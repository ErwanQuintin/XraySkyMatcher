# XraySkyMatcher
A light version of the STONKS package, meant to allow to do cone-searches in the X-ray multi-instrument catalog, and retrieve the corresponding sources' properties

Requirements:
- STILTS (https://www.star.bris.ac.uk/~mbt/stilts/)
- Standard Python packages

"MatchCatalogs.py" will use STILTS to do a cone search around the target position with the input radius. It will also crop the instrument catalogs to only contain the relevant sources.
"LoadMasterSources.py" will load the relevant multi-instrument sources in custom Python objects. This gives access to flux and band photometry.

To function, the catalog data needs to be downloaded from the following link: https://drive.google.com/file/d/1ZdYlraZQXDN_grrz4aGwOFBCO_8nKidZ/view?usp=drive_link and unzipped in the same file as the Python scripts.
