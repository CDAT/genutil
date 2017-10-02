"""This file contains docstring snippets that are used repetetively in documenting genutil functions.
    Use it as a repository for any pieces of documentation that show up many times, and would be more convenient to
    edit in one place and plug in everywhere.
"""

# Parameters
axis = """
    :param axis: 'x' | 'y' | 'z' | 't' | '(dimension_name)' | 0 | 1 ... | n
        default value = 0. You can pass the name of the dimension or index
        (integer value 0...n) over which you want to compute the statistic.
        you can also pass 'xy' to work on both axes at once
    :type axis: str or int
    """

df = """
    :param df: An integer flag indicating whether to return degrees of freedom.
        If 1, degrees of freedom are returned. If set to 0, they are not.
    :type df: int
    """
