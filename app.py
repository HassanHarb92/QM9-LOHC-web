from flask import Flask, render_template, request
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Draw
import numpy as np
import os
from io import BytesIO
import base64

app = Flask(__name__)


# Load the dataset
#QM9_G4MP2_all = pd.read_csv('QM9_G4MP2_all.csv')
QM9_G4MP2_all = pd.read_csv('G4MP2_set2.csv')

@app.route('/', methods=['GET', 'POST'])
def index():
    if request.method == 'POST':
        # Extract parameters from form data
        print(request.form)
        delta_H_range = (float(request.form['delta_H_min']), float(request.form['delta_H_max']))
        pH2_range = (float(request.form['pH2_min']), float(request.form['pH2_max']))
        molar_mass_range = (0.0,9999.99) # (float(request.form['molar_mass_min']), float(request.form['molar_mass_max']))
        num_results = int(request.form['num_results'])
        
        # Filter the dataset based on the specified ranges
        filtered_data = QM9_G4MP2_all[
            (QM9_G4MP2_all['delta_H'].between(*delta_H_range)) &
            (QM9_G4MP2_all['pH2'].between(*pH2_range)) #&
#            (QM9_G4MP2_all['molar_mass'].between(*molar_mass_range))
        ]
        
#        if len(filtered_data) > 10:
#            filtered_data = filtered_data.sample(n=3)

        # Limit the number of results if there are more results than the user requested
        if len(filtered_data) > num_results:
            filtered_data = filtered_data.sample(n=num_results)
        
        # Convert molecule images to base64 for embedding in HTML
        def mol_to_image_base64(smiles):
            mol = Chem.MolFromSmiles(smiles)
            img = Draw.MolToImage(mol)
            output = BytesIO()
            img.save(output, format='PNG')
            return base64.b64encode(output.getvalue()).decode('utf-8')
        
        filtered_data['molecule1_img'] = filtered_data['unsat_SMILE'].apply(mol_to_image_base64)
        filtered_data['molecule2_img'] = filtered_data['sat_SMILE'].apply(mol_to_image_base64)
        actual_num_results = len(filtered_data)  # Get the actual number of results after filtering and possibly sampling

        return render_template('results.html', data=filtered_data, num_results=actual_num_results)
        

        # Pass the filtered data to your template
#        return render_template('results.html', tables=[filtered_data.to_html(classes='data')], titles=filtered_data.columns.values, data=filtered_data)
    
    # If method is GET, just render the form
    return render_template('index.html')



