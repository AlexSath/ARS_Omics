def n_underscore_constructor(n, keyword_dictionary):
    def this_n_underscore(path):
        spot = path.split("_")[n]
        try:
            return keyword_dictionary[spot]
        except KeyError:
            return spot
    return this_n_underscore