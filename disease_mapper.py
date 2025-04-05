import requests
import time
import os


def map_diseases(disease_terms, api_key=None):
    if not api_key:
        api_key = os.environ.get("BIOONTOLOGY_API_KEY", "")

    if not api_key or api_key == "":
        raise ValueError(
            "API key is required. Set the BIOONTOLOGY_API_KEY environment variable or pass it as an argument."
        )

    mondo_dict = {}
    doid_dict = {}
    mesh_dict = {}

    ontologies = ["MONDO", "DOID", "MESH"]
    ontologies_str = ",".join(ontologies)

    for term in disease_terms:
        print(f"Processing: {term}")
        url = "https://data.bioontology.org/search"

        params = {
            "q": term,
            "ontologies": ontologies_str,
            "require_exact_match": "false",
            "include": "prefLabel",
            "apikey": api_key,
        }

        response = requests.get(url, params=params)

        if response.status_code == 200:
            data = response.json()

            best_matches = {
                "MONDO": {"score": 0, "id": None},
                "DOID": {"score": 0, "id": None},
                "MESH": {"score": 0, "id": None},
            }

            if data["collection"]:
                for item in data["collection"]:
                    ontology = item["links"]["ontology"].split("/")[-1]

                    if ontology not in ontologies:
                        continue

                    full_id = item["@id"]
                    local_id = full_id.split("/")[-1]

                    pref_label = item.get("prefLabel", "")

                    if term.lower() == pref_label.lower():
                        match_score = 1.0

                    elif (
                        term.lower() in pref_label.lower()
                        or pref_label.lower() in term.lower()
                    ):
                        match_score = 0.7
                    else:
                        match_score = 0.3

                    if match_score > best_matches[ontology]["score"]:
                        best_matches[ontology]["score"] = match_score
                        best_matches[ontology]["id"] = local_id

            if best_matches["MONDO"]["id"]:
                mondo_dict[term] = best_matches["MONDO"]["id"]
            if best_matches["DOID"]["id"]:
                doid_dict[term] = best_matches["DOID"]["id"]
            if best_matches["MESH"]["id"]:
                mesh_dict[term] = best_matches["MESH"]["id"]

        else:
            print(f"Error for term '{term}': {response.status_code}")

        time.sleep(0.5)

    return mondo_dict, doid_dict, mesh_dict
