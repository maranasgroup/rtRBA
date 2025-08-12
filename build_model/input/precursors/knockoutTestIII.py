import pandas as pd
import cobra
from collections import OrderedDict
import copy
from copy import deepcopy
import os
from itertools import combinations

curr_dir = os.getcwd()

# Metabolic model (COBRApy json)
model = cobra.io.load_json_model(curr_dir + '/../iRhtoCN.json')
biomass = 'BIOMASS'

# Define the set of reactions to knock out
#all_reactions = set(['CHTNDA_c','RPE_c','SDPDS_c','TKT1_c','TKT2_c','XU5PFGT_x'])
all_reactions = set(['CHTNDA_c','RPE_c','TKT2_c','XU5PFGT_x'])

# Function to perform FBA and categorize reactions
def categorize_reactions(model, reactions_to_knockout, essential_reactions, feasible_knockouts):
    # Create a copy of the model to avoid modifying the original model
    model_copy = model.copy()

    # Knock out the reactions in the current combination
    for reaction_id in reactions_to_knockout:
        reaction = model_copy.reactions.get_by_id(reaction_id)
        reaction.lower_bound = 0.0
        reaction.upper_bound = 0.0

    # Perform FBA on the modified model
    model_copy.objective = model.reactions.get_by_id(biomass)
    model_copy.objective_direction = "max"
    knockout_solution = model_copy.optimize()

    # Check the objective value
    if knockout_solution.objective_value <= 0:
        essential_reactions.update(reactions_to_knockout)
        print("Lethal knockout: ", reactions_to_knockout)
    else:
        feasible_knockouts.extend(reactions_to_knockout)
        print("Nonlethal knockout: ", reactions_to_knockout)

# Lists to store results
essential_reactions = set()
feasible_single_knockouts = []
feasible_knockout_combinations = [feasible_single_knockouts]

# Try single knockouts
for reaction_id in all_reactions:
    categorize_reactions(model, {reaction_id}, essential_reactions, feasible_single_knockouts)

# Function to generate higher-order knockout combinations
def generate_higher_order_combinations(feasible_knockouts, n):
    higher_order_combinations = []
    for combo in combinations(feasible_knockouts, n):
        # Sort reactions within the combination to ensure consistency
        combined_knockout = frozenset(sorted(item for sublist in combo for item in sublist))
        higher_order_combinations.append(combined_knockout)
    return higher_order_combinations

# Try higher-order knockouts (e.g., double, triple, quadruple)
for order in range(2, len(all_reactions)):
    higher_order_combinations = generate_higher_order_combinations(feasible_knockout_combinations[-1], order)
    feasible_knockout_combinations.append(set())  # Use sets to ensure uniqueness
    for combo in higher_order_combinations:
        categorize_reactions(model, combo, essential_reactions, feasible_knockout_combinations[-1])

# Print the results
print("Essential Reactions:")
print(essential_reactions)
print("\nFeasible Single Knockouts:")
print(feasible_single_knockouts)
print("\nFeasible Knockout Combinations:")
for i, combinations in enumerate(feasible_knockout_combinations):
    print(f"Order-{i + 2} Knockouts:")
    print(combinations)
