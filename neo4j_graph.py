import pandas as pd
from neo4j import GraphDatabase
import numpy as np
import dotenv
import os


load_status = dotenv.load_dotenv("data/Neo4j-e6fae6f2-Created-2025-04-05.txt")
if load_status is False:
    raise RuntimeError('Environment variables not loaded.')

URI = os.getenv("NEO4J_URI")
AUTH = (os.getenv("NEO4J_USERNAME"), os.getenv("NEO4J_PASSWORD"))

with GraphDatabase.driver(URI, auth=AUTH) as driver:
    driver.verify_connectivity()
    print("Connection established.")