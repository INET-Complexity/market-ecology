import esl

print(esl.version())




class SingleExperiment(esl.simulation.model):
    def __init__(self, env, par):
        super().__init__(self, env,par)




if __name__ == "__main__":
    e = esl.computation.environment()
    esl.simulation.parameter.par
    m = SingleExperiment(e, {})




def clear_market(initial_prices, excess_demand_functions):
    model = excess_demand_model(initial_prices)
    model.excess_demand_functions = excess_demand_functions
    return model.compute_clearing_quotes()
