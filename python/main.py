import esl

class SingleExperiment(esl.simulation.model):
    def __init__(self, environment, parameters):
        super().__init__(environment, parameters)

        self.economy = esl.simulation.world()




if __name__ == "__main__":
    e = esl.computation.environment()
    ps = esl.simulation.parameter.parametrization()
    ks = ps.values.keys()
    print(ks[0])
    m = SingleExperiment(e, ps)

    e.run(m)

def clear_market(initial_prices, excess_demand_functions):
    model = excess_demand_model(initial_prices)
    model.excess_demand_functions = excess_demand_functions
    return model.compute_clearing_quotes()
