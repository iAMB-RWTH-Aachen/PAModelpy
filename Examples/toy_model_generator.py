from Scripts.toy_model_TS_generation import build_toy_model

toy_mcpam = build_toy_model(sensitivity=False)
toy_mcpam.objective = 'R11'
result = toy_mcpam.optimize()
print(toy_mcpam.objective.value)
