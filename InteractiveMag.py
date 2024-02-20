import matplotlib.pyplot as plt
from planetmagfields import Planet, get_models, utils
import re

# Assuming get_models and Planet are defined elsewhere in your code or imported from a library

# def extract_planet_names_from_file() -> list:
#     """
#     Extracts the names of planets from specified .rst file within the local 'doc' directory.

#     This function reads the 'models.rst' file, searching for planet names using a regular expression pattern.
#     The names are expected to be in a specific format within the file, marked by asterisks and colons.

#     Returns:
#         list of str: A list containing the names of the planets extracted from the file.
#     """
#     file_path = "doc/models.rst"  # Relative path from the notebook to the models.rst file

#     # Read the content of the models.rst file
#     with open(file_path, 'r') as file:
#         content = file.read()

#     # Regex pattern to accurately capture planet names
#     pattern = r"\*\s*\*(.*?)\*\s*:"

#     # Find all matches in the content
#     planet_matches = re.findall(pattern, content)

#     # List to hold the names of the planets
#     planets = [planet.strip() for planet in planet_matches]

#     return planets

def extract_models_for_planets() -> list:
    """
    Extracts and modifies model names for a given list of planet names.

    For each planet name provided, this function retrieves associated model names,
    modifies these model names by prepending the planet's name, and then aggregates
    all modified model names into a single list.

    Args:
        planet_names (list of str): A list of planet names for which to extract and modify model names.

    Returns:
        list of str: A list containing all modified model names for the given planet names.
    """
    all_model_names = []  # Initialize an empty list to store all model names

    for planet_name in utils.planetlist:
        model_names = get_models(planet_name)

        # Prepend planet name to each model name
        # modified_model_names = [f"{planet_name}_{model_name}" for model_name in raw_model_names]

        # Print all modified models for reference
        print(planet_name, ":", model_names)

        # Extend the all_model_names list with the modified model names
        all_model_names.extend(model_names)

    return all_model_names

# def extract_only_model_names(array: list) -> list:
#     """
#     Extracts only the model names from a list of strings formatted as "planetName_modelName".

#     This function processes each string in the input list, removing the planet name and underscore,
#     leaving only the model name. If a model name contains underscores, they are preserved.

#     Args:
#         array (list of str): A list of strings, each formatted as "planetName_modelName".

#     Returns:
#         list of str: A list of model names with the planet names removed.
#     """
#     model_names = ['_'.join(item.split('_')[1:]) for item in array]
#     return model_names

def plot_intercat_mag_r(name: str, r: float, model: str, background: str) -> None:
    """
    Plots the magnetic field of a planet for a given model and radius, with a specified background.

    This function sets the plotting style based on the background preference ('Dark' or 'Default'),
    creates a planet object with the specified name and model, and then plots the planet's magnetic
    field at a given radius in a Mollweide projection. It also displays the spectral characteristics
    of the magnetic field.

    Args:
        name (str): The name of the planet.
        r (float): The radial distance at which to plot the magnetic field.
        model (str): The model name to use for plotting.
        background (str): The background color style for the plot ('Dark' or 'Default').

    Returns:
        Intervactive Plot
    """
    # Set background color based on user choice
    if background == 'Dark':
        plt.style.use('dark_background')
    else:
        plt.style.use('default')

    p = Planet(name=name, model=model, datDir=utils.stdDatDir)
    p.plot(r=r)
    plt.show()

    p.spec(r=r)
    plt.show()

def save_figure_with_options(title: str = 'figure', dpi: int = 600, transparent: bool = False) -> None:
    """
    Saves the current matplotlib figure to a PNG file with customizable options for the title, resolution (DPI), and background transparency.

    This function is designed to be used in a Jupyter notebook or any Python environment where matplotlib figures are created. It allows the user to easily save the current figure with specified options directly from their code, enhancing the workflow for generating and exporting figures for reports, presentations, or further analysis.

    Parameters:
        title (str): The title and filename under which to save the figure. The default filename is 'figure'.
        dpi (int): The resolution of the saved figure in dots per inch (DPI). Higher DPI values result in higher resolution images. The default DPI is 600.
        transparent (bool): A flag indicating whether the figure should be saved with a transparent background. If True, the figure's background will be transparent. Otherwise, it will have the default opaque background. The default is False.

    Returns:
        None: This function does not return a value but saves the current figure to a PNG file with the specified options.

    Examples:
        # To save the figure with default settings (title='figure', dpi=300, transparent=False):
        save_figure_with_options()

        # To save the figure with a custom title, increased DPI, and transparent background:
        save_figure_with_options(title='MyCustomFigure', dpi=600, transparent=True)
    """
    plt.savefig(f'{title}.png', dpi=dpi, transparent=transparent)
    print(f"Figure saved as '{title}.png' with dpi={dpi} and transparent={transparent}")

