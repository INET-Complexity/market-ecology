import esl

class SingleExperiment(esl.model):
    def __init__(self):
        # SingleExperiment()
        pass



if __name__ == "__main__":

    print(esl.version())




def clear_market(initial_prices, excess_demand_functions):
    model = excess_demand_model(initial_prices)
    model.excess_demand_functions = excess_demand_functions
    return model.compute_clearing_quotes()
