
import streamlit as st
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Draw
import numpy as np
import base64
from PIL import Image
import io

# Load the dataset
# QM9_G4MP2_all = pd.read_csv('QM9_G4MP2_all.csv')
QM9_G4MP2_all = pd.read_csv('G4MP2_set2.csv')

# Convert molecule images to base64 for embedding
def mol_to_image_base64(smiles):
    mol = Chem.MolFromSmiles(smiles)
    img = Draw.MolToImage(mol)
    output = io.BytesIO()
    img.save(output, format='PNG')
    return base64.b64encode(output.getvalue()).decode('utf-8')

# Streamlit app
st.title('My Flask App Converted to Streamlit')

delta_H_min = st.number_input('Delta H Min', value=float(QM9_G4MP2_all['delta_H'].min()))
delta_H_max = st.number_input('Delta H Max', value=float(QM9_G4MP2_all['delta_H'].max()))
pH2_min = st.number_input('pH2 Min', value=float(QM9_G4MP2_all['pH2'].min()))
pH2_max = st.number_input('pH2 Max', value=float(QM9_G4MP2_all['pH2'].max()))
num_results = st.number_input('Number of Results', value=10, min_value=1, max_value=len(QM9_G4MP2_all))

if st.button('Submit'):
    filtered_data = QM9_G4MP2_all[
        (QM9_G4MP2_all['delta_H'].between(delta_H_min, delta_H_max)) &
        (QM9_G4MP2_all['pH2'].between(pH2_min, pH2_max))
    ]

    if len(filtered_data) > num_results:
        filtered_data = filtered_data.sample(n=num_results)

    filtered_data['molecule1_img'] = filtered_data['unsat_SMILE'].apply(mol_to_image_base64)
    filtered_data['molecule2_img'] = filtered_data['sat_SMILE'].apply(mol_to_image_base64)

    st.write(filtered_data[['unsat_SMILE', 'sat_SMILE', 'molecule1_img', 'molecule2_img']])

    for index, row in filtered_data.iterrows():
        col1, col2 = st.columns(2)
        with col1:
            st.image(Image.open(io.BytesIO(base64.b64decode(row['molecule1_img']))), caption="Unsaturated")
        with col2:
            st.image(Image.open(io.BytesIO(base64.b64decode(row['molecule2_img']))), caption="Saturated")

