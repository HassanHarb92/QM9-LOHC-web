import streamlit as st
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Draw
from PIL import Image
import io
import base64

# Assuming 'G4MP2_set2.csv' is present in the same directory as your Streamlit app
QM9_G4MP2_all = pd.read_csv('G4MP2_set2.csv')

# Convert molecule images to base64 for embedding
def mol_to_image_base64(smiles):
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is not None:
            img = Draw.MolToImage(mol)
            output = io.BytesIO()
            img.save(output, format='PNG')
            return base64.b64encode(output.getvalue()).decode('utf-8')
        else:
            return None
    except Exception as e:
        st.error(f"Error generating image for SMILES: {smiles} | Error: {str(e)}")
        return None

#st.title('My Flask App Converted to Streamlit')
st.title('QM9-LOHC Dataset Query')

image_path = 'Logo.png'
st.image(image_path, caption="hello", use_column_width=True)



# User inputs for filtering
delta_H_min = st.number_input('Delta H Min', value=0)# float(QM9_G4MP2_all['delta_H'].min()))
delta_H_max = st.number_input('Delta H Max', value=0)#float(QM9_G4MP2_all['delta_H'].max()))
pH2_min = st.number_input('%H2 (by wt) Min', value=0)#float(QM9_G4MP2_all['pH2'].min()))
pH2_max = st.number_input('%H2 (by wt) Max', value=0)#float(QM9_G4MP2_all['pH2'].max()))
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

    st.write(filtered_data[['unsat_SMILE', 'sat_SMILE','pH2']])#, 'molecule1_img', 'molecule2_img']])

    # Displaying images
    for index, row in filtered_data.iterrows():
        col1, col2 = st.columns(2)
        with col1:
            if row['molecule1_img'] is not None:
                st.image(Image.open(io.BytesIO(base64.b64decode(row['molecule1_img']))), caption="Unsaturated")
        with col2:
            if row['molecule2_img'] is not None:
                st.image(Image.open(io.BytesIO(base64.b64decode(row['molecule2_img']))), caption="Saturated")

