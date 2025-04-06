import requests
import pandas as pd
from tqdm import tqdm
import time

# NOTE: This code was adapted from https://clarewest.github.io/blog/post/crash-course-in-open-targets-part-1/
class OpenTargetsClient:
    
    def __init__(self, base_url="https://api.platform.opentargets.org/api/v4/graphql"):
        self.base_url = base_url
        
    def execute_query(self, query, variables=None):
        payload = {"query": query}
        if variables:
            payload["variables"] = variables
            
        response = requests.post(self.base_url, json=payload)
        
        if response.status_code != 200:
            raise Exception(f"Query failed with status code {response.status_code}: {response.text}")
            
        return response.json()
    
    def get_disease_associated_targets(self, efo_id, limit=100):
        query = """
        query DiseaseAssociatedTargets($efoId: String!, $limit: Int!) {
          disease(efoId: $efoId) {
            id
            name
            associatedTargets(page: {index: 0, size: $limit}) {
              count
              rows {
                target {
                  id
                  approvedSymbol
                  approvedName                
                }
                score
                datatypeScores {
                  id
                  score
                }
              }
            }
          }
        }
        """
        
        variables = {
            "efoId": efo_id,
            "limit": limit
        }
        
        return self.execute_query(query, variables)


def get_targets_for_efo_ids(efo_ids, output_file=None, rate_limit=1):
    client = OpenTargetsClient()
    all_results = []
    
    for efo_id in tqdm(efo_ids, desc="Processing EFO IDs"):
        try:
            response = client.get_disease_associated_targets(efo_id, limit=500)
            
            disease_data = response.get('data', {}).get('disease', {})
            
            if not disease_data:
                print(f"No data found for EFO ID: {efo_id}")
                continue
                
            disease_name = disease_data.get('name', 'Unknown')
            disease_id = disease_data.get('id', efo_id)
            
            target_rows = disease_data.get('associatedTargets', {}).get('rows', [])
            
            if not target_rows:
                print(f"No target associations found for {disease_name} ({efo_id})")
                continue
                
            print(f"Found {len(target_rows)} target associations for {disease_name} ({efo_id})")
            
            for row in target_rows:
                target = row.get('target', {})
                
                result = {
                    'disease_id': disease_id,
                    'disease_name': disease_name,
                    'ensembl_id': target.get('id'),
                    'approved_symbol': target.get('approvedSymbol'),
                    'approved_name': target.get('approvedName'),
                    'association_score': row.get('score')
                }
                
                all_results.append(result)
            
            time.sleep(rate_limit)
            
        except Exception as e:
            print(f"Error processing EFO ID '{efo_id}': {str(e)}")
            time.sleep(rate_limit)

    if all_results:
        result_df = pd.DataFrame(all_results)
        
        if output_file:
            result_df.to_parquet(output_file, index=False)
            print(f"Saved {len(result_df)} associations to {output_file}")
            
        return result_df
    else:
        print("No associations found")
        return None
