from flask import Flask,request,jsonify
from flask_cors import CORS
import os
from dotenv import load_dotenv
import main

load_dotenv()

app = Flask(__name__)

# Access the environment variables
app.config['FRONTEND_URL'] = os.getenv('FRONTEND_URL')
# Setup CORS
CORS(app, origins=["http://localhost:5173",app.config['FRONTEND_URL']])

# simulation route
@app.route('/simulate', methods=['POST'])
def get_transient_res():

    data = request.json
    circuit = data.get('circuit')
    simulationType = data.get('simulationType')
    simulationTimeStep = data.get('simulationTimeStep') or 0.01
    stamps = data.get('stamps') or 10

    result = main.main(circuit,simulationType,simulationTimeStep,stamps)

    return jsonify(message="All is well", result=result), 200
    
# test route
@app.route('/test', methods=['GET'])
def get_items():
    return jsonify({"message": "success"}), 201

# Run the app
if __name__ == '__main__':
    app.run(debug=True)