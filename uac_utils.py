# Import necessary modules
import re
import importlib
import matplotlib.pyplot as plt
import pandas as pd
import geopandas as gpd
from shapely.geometry import box, Polygon, MultiPolygon, LineString
from shapely.geometry import MultiLineString
from shapely.ops import polygonize, unary_union
from pyproj import CRS
import rasterio

def reproject_gdf(gdf, epsg_code):
    """Reprojects GeoDataFrame to CRS with EPSG code

    Assigns WKT format of projection to EPSG code to GeoDataFrame.

    Args:
    gdf: GeoDataFrame with any geometry (e.g., Point, Line, Polygon)
    epsg_code: EPSG code (integer)

    Returns:
        GeoDataFrame reprojected to new crs (based on EPSG code).
    """
    # Define CRS in WKT format using EPSG code
    target_projection = CRS.from_epsg(epsg_code).to_wkt()

    # Reproject GeoDataFrame to epsg_code
    reprojected_gdf = gdf.to_crs(target_projection)

    # Print message
    print("GeoDataFrame now has the following CRS:\n")
    print(reprojected_gdf.crs)

    return reprojected_gdf

def gdf_has_duplicate_rows(gdf):
    """Returns True if gdf GeoDataFrame has duplicate rows

    Args:
        gdf: GeoDataFrame

    Returns:
        Boolean value returned based on whether gdf has duplicate
        rows (True) or not (False)

    """

    # Create mask for all duplicate rows
    gdf_duplicate_mask = gdf.duplicated()

    # Calculate # duplicate rows
    number_duplicate_rows = len(gdf[gdf_duplicate_mask])

    # Return True if there are any duplicate rows
    return number_duplicate_rows > 0

def create_index_column(gdf, index_colname='index'):
    """Returns GeoDataFrame with index as separate column

    Args:
        gdf: GeoDataFrame
        index_colname: name of new column that holds index.
            Set to 'index' by default

    Returns:
        gdf_with_index, a GeoDataFrame with index as separate column.
    """

    # Make copy of gdf
    gdf_with_index = gdf.copy()

    # Create new column with index
    gdf_with_index[index_colname] = gdf_with_index.index

    return gdf_with_index

def all_polygon_geometries(gdf):
    """Returns True if all geometries are (Multi)Polygon"""

    # Make copy of gdf
    gdf_copy = gdf.copy()

    # Create new column with type of each row geometry
    gdf_copy['geom_type'] = type(gdf_copy['geometry'])

    # Find unique values of geometry types
    geom_type_list = gdf_copy.geom_type.unique()

    # Check if each type is a Polygon feature
    geom_is_poly = ["Polygon" in geom for geom in geom_type_list]

    # If at least one geometry is non-Polygon, return False
    # Otherwise, return True
    if False in geom_is_poly:
        return False
    else:
        return True

def create_map_registration_dict(gdf, map_colname='MAP_NO',
                                 registration_colname='REGISTRATI'):
    """Create dictionary with map number (key) and registration number (value)

    Args:
        gdf: GeoDataFrame with Map
        map_colname: Name of column with map no. Since this name varies
            slightly with different shapefiles, it is an optional
            parameter, set to 'MAP_NO' as default.
        registration_colname: Name of column with registration no.
            Since this name varies slightly with different shapefiles,
            it is an optional parameter, set to 'MAP_NO' as default.

    Returns:
        A dictionary with the map number as the key and the
        registration number as the value.
    """
    # Generate unique map number
    unique_map_numbers = gdf[map_colname].unique()

    # Initialize empty dictionary
    map_registration_dict = dict()

    # Iterate across all unique map numbers
    for map_number in unique_map_numbers:
        # Create a set object
        map_registration_dict[map_number] = set()

        # Iterate over all rows where Map Number equals map_number
        for idx, row in gdf[gdf[map_colname] == map_number].iterrows():

            # Add Registration number to set
            map_registration_dict[map_number].add(row[registration_colname])

    return map_registration_dict

def maps_numbers_with_duplicate_polygons(gdf, map_colname='MAP_NO'):
    """Create list of map numbers that are duplicated across polygons

    Args:
        gdf: GeoDataFrame with Map
        map_colname: Name of column with map no. Since this name varies
            slightly with different shapefiles, it is an optional
            parameter, set to 'MAP_NO' as default.
    Returns:
        A list of map numbers that are assigned to multiple polygons.
    """
    duplicate_mask = gdf[map_colname].duplicated()

    duplicate_map_numbers = gdf[duplicate_mask][map_colname].unique()

    print("Number of map numbers that are duplicated: {}".\
        format(len(duplicate_map_numbers)))

    return duplicate_map_numbers

def create_map_reg_codes(map_registration_dict, duplicate_map_numbers):
    """Create dictionary that associates map number with string for PDF URLs

    This function creates a dictionary that has the following structure for
    each map number: {map_number: '{map_number}_{registration_number}'}. This
    string is a distinguishing feature of every URL for a georeferenced PDF
    map.

    Args:
        map_registration_dict: Dictionary that associates each map number
            with registration number.
        duplicate_map_numbers: List of map numbers associated with more
            than one polygon

    Returns:
        Dictionary of the form: {map_number:
        '{map_number}_{registration_number}'}
    """

    # Create empty dictionary
    map_reg_code_dict = dict()

    # Form new string: '{map_number}_{registration_number}' and associate this
    # with each map number
    map_reg_code_dict = {map_num:"{}_{}".format(map_num,
                            list(map_registration_dict[map_num])[0])\
                            for map_num in duplicate_map_numbers}

    return map_reg_code_dict

def clean_map_reg_string(map_reg_string):
    """Remove trailing characters from mapreg_string

    Removes extraneous whitespace and other characters that
    are part of the Map No/Registration No string that are
    not likely part of the PDF URL.

    Args:
        map_reg_string: String that concatenates map number
            with registration number.

    Returns:
        String that removes whitespace, dash (-), and opening
        parenthesis that follow the "mapno_regno" part
        of the string.

    """
    # Find first occurrence of single space, dash, and
    # open parenthesis. If any of these is not found,
    # the .find() method returns -1.
    index1 = map_reg_string.find(' ')
    index2 = map_reg_string.find('-')
    index3 = map_reg_string.find('(')

    # Initialize clean_map_reg_string as input string.
    # In case none of the characters are found, we return
    # the original string input.
    clean_map_reg_string = map_reg_string

    # Create list of indices
    idx_list = [index1, index2, index3]

    # Find all non-negative numbers
    non_neg_indices = [num for num in idx_list if num>=0]

    # If there is non-negative index, find the minimum index
    # and use this to subset string
    if non_neg_indices:
        return clean_map_reg_string[:min(non_neg_indices)]
    # otherwise, return the input string
    else:
        return clean_map_reg_string

def extract_map_URL_segment(map_number, dda_file):
    """
    Extracts URL segment for Georeferenced PDF of a specific map number

    Args:
        map_number: Map number for which we want PDF URL.
        dda_file: String with the DDA html content

    Returns:
        String with the format:
        '/tendernotices_docs/27112019/{mapno}_{regno}.pdf' OR
        '/tendernotices_docs/27112019/MAP {mapno}_{regno}.pdf'
        where mapno is the map number and regno is the registration
        number. Note that some URLs will add characters after the
        registration number and that regex pattern accounts for this.
    """

    # Define pattern to search for map number PDF URL
    # Parentheses allow extraction of group/string we want
    pattern = r"\'..(/tendernotices_docs/\d+/{}_.*.pdf)".format(map_number)

    # Do regex matching on DDA HTML content
    matches = re.search(pattern, dda_file)

    url_prefix = r"https://dda.org.in/"

    # If a match is found, return URL path
    # Replace any whitespace character with "%20"
    if matches:
        return url_prefix + matches.group(1).replace(" ", "%20")
    else:
        # Try another pattern for the URL
        pattern = r"..(/tendernotices_docs/\d+/MAP {}_.*.pdf)".\
            format(map_number)
        matches = re.search(pattern, dda_file)

        if matches:
            # Return URL path, replacing whitespace with "%20"
            return url_prefix + matches.group(1).replace(" ", "%20")
        else:
            return "URL not found"

def generate_map_URLs(map_number_list):
    """Generate URLs from map numbers of unauthorized colonies

    Args:
        map_number_list: list of map numbers for which there are
            duplicate polygons associated with them

    Returns:
        Dictionary of the form: {map_number: map_pdf_url},
        where map_number is an integer and map_pdf_url
        is similar to the following string:
        "https://dda.org.in/tendernotices_docs/27112019/733_291.pdf" or
        "https://dda.org.in//tendernotices_docs/27112019/MAP%20202_594.pdf"
    """

    # Extract DDA HTML content to variable dda_file
    dda_html_filepath = "DDA_webpage_html.txt"
    with open(dda_html_filepath) as f:
        dda_file = f.read()

    # Initialize empty dictionary
    map_url_dict = dict()

    # Extract URL for each map number and add the map number
    # and map URL to map_url_dict dictionary.
    for map_number in map_number_list:
        map_url = extract_map_URL_segment(map_number, dda_file)
        map_url_dict[map_number] = map_url

    return map_url_dict

def generate_map_centroid_crs(map_url_dict):
    """Create dictionary associating map number with centroid and CRS

    Args:
        map_url_dict: Dictionary with map number as key and map
            PDF url as value.

    Returns:
        Dictionary with the following format for each entry:
        {map number: {'centroid': POINT(...), 'crs': 'EPSG:ABCDE'}}.
    """

    # Initialize map_centroid_crs_dict
    map_centroid_crs_dict = dict()

    # Iterate over elements of map_url_dict
    for map_number, map_url in map_url_dict.items():
        # Extract bounding box and CRS from Georeferenced PDF
        with rasterio.open(map_url) as geopdf:
            pdf_bounds = geopdf.bounds
            pdf_crs = geopdf.crs

        # Compute centroid of PDF
        pdf_centroid = box(pdf_bounds[0], pdf_bounds[1],
                           pdf_bounds[2], pdf_bounds[3]).centroid

        # Add map number, centroid and crs to map_centroid_crs_dict
        map_centroid_crs_dict[map_number] = {"centroid": pdf_centroid,
                                            "crs": pdf_crs}

    return map_centroid_crs_dict

def identify_polygon_from_map_no_and_pdf(uac_gdf, map_number,
                                            map_centroid_crs_dict,
                                            map_colname='MAP_NO'):
    """Identifies which polygon is most likely associated with map_no

    Finds all polygons associated with a specific map_number and computes
    distances from their centroids to the centroid of the bounding box of the
    PDF from which these polygons were extracted. The polygon with the shortest
    distance to the bounding box centroid is selected.

    Args:
        uac_gdf: GeoDataFrame with polygon geometries. It is assumed
            to be the data for unauthorized colonies.
        map_number: Map_No for which there are duplicate entries in
            polygon_gdf.
        map_centroid_crs_dict: Dictionary with the following
            format for each entry:
            {map number: {'centroid': POINT(...), 'crs': 'EPSG:ABCDE'}}.
        map_colname: Name of column with map no. Since this name varies
            slightly with different shapefiles, it is an optional
            parameter, set to 'MAP_NO' as default.

    Returns:
        polygons_to_keep_delete: Python distionary with two keys: "keep"
            and "delete". "keep" has as its value the index of the
            polygon to keep while "delete" has as its value a list
            of indices for all rows to delete.
    """
    # Create a subset of uac_gdf only consisting of rows
    # where "Map_No" equals map_number
    uac_gdf_subset =  uac_gdf[uac_gdf[map_colname] == map_number]

    # Retrieve centroid and CRS of PDF
    pdf_crs = map_centroid_crs_dict[map_number]['crs']
    pdf_centroid = map_centroid_crs_dict[map_number]['centroid']

    # Reproject uac_gdf_subset to PDF's CRS
    uac_gdf_subset =  uac_gdf_subset.to_crs(pdf_crs)

    # Initialize dictionary with centroid distances
    # Key is index and value is distance
    centroid_distances = dict()

    # Iterate across all rows in subset matching map_number
    for idx, row in uac_gdf_subset.iterrows():
        # Add centroid distances to centroid_distance dictionary
        centroid_distances[row['index']] = row['geometry'].centroid.\
                                            distance(pdf_centroid)

    # Calculate the minimum distance of all centroid distances
    min_distance = min(centroid_distances.values())

    # Initialize dictionary with polygons to keep and delete
    polygons_to_keep_delete = dict()

    # Initialize delete key with a value of an empty list
    polygons_to_keep_delete['delete'] = []

    # If index corresponds to shortest distance designate as "keep"
    # else place the index list for the "delete" key
    for index in centroid_distances.keys():
        if centroid_distances[index] == min_distance:
            polygons_to_keep_delete['keep'] = index
        else:
            polygons_to_keep_delete['delete'].append(index)

    return polygons_to_keep_delete

def nearest_polygon_border(uac_gdf, map_number, map_centroid_crs_dict,
                            map_colname="MAP_NO"):
    """Find nearest polygon by distance from PDF centroid to border"""

    # Create a subset of uac_gdf only consisting of rows
    # where "Map_No" equals map_number
    uac_gdf_subset = uac_gdf[uac_gdf[map_colname] == map_number]

    # Retrieve centroid and CRS of PDF
    pdf_crs = map_centroid_crs_dict[map_number]['crs']
    pdf_centroid = map_centroid_crs_dict[map_number]['centroid']

    # Reproject uac_gdf_subset to PDF's CRS
    uac_gdf_subset =  uac_gdf_subset.to_crs(pdf_crs)

    nearest_distance_dict = dict()

    for idx, row in uac_gdf_subset.iterrows():
        nearest_polygon_point = nearest_points(pdf_centroid, row['geometry'])[1]
        nearest_distance = nearest_polygon_point.distance(pdf_centroid)
        polygon_centroid_distance = row['geometry'].centroid.\
                                    distance(pdf_centroid)
        nearest_distance_dict[row['index']] = {'nearest_distance':\
                                                nearest_distance,
                                              'polygon_centroid_distance':\
                                               polygon_centroid_distance}

    #nearest_indices = [key for key, val in nearest_distance_dict.items() \
    #                   if val==min(nearest_distance_dict.values())]

    #if len(nearest_indices) > 1:
    #    return "More than 1"
    #else:
    #    return nearest_indices

    return nearest_distance_dict

## Convert LineStrings to Polygons (Disbanded)
#    * Check if geometry is not (Multi)Polygon
#    * Apply `shapely.ops.polygonize` to (Multi)LineString
#    * Merge all (Multi)Polygon objects as one merged polygon
#    * Return merged (Multi)Polygon object

def convert_linestring_to_polygon(geom):

    """Convert (Multi)LineString to (Multi)Polygon

    Use shapely.ops.polygonize to transform (Multi)LineString
    into Multi(Polygon). The function merges all (Multi)Polygons
    into one geometry.

    Args:
        geom: Multi(LineString) geometry

    Returns:
        merged_geom: Multi(Polygon) geometry
    """

    # Check that geometry is (Multi)LineString
    assert isinstance(geom, (LineString, MultiLineString)),\
        "geometry is not (Multi)LineString"

    # Use shapely.ops.polygonize to transform
    # (Multi)LineString into (Multi)Polygon(s)
    # Note that the output is a generator yielding
    # (Multi)Polygon geometries
    polygonized_geom = polygonize(geom)

    # Merge the Multi(Polygon) geometry(ies)
    merged_geom = unary_union(list(polygonized_geom))

    return merged_geom

def convert_geoms_to_polygon(gdf):
    """Check/convert all gdf geometries to (Multi)Polygon

    Checks all geometries in a GeoDataFrame gdf that should consist
    of (Multi)Polygons but may have (Multi)LineStrings. If the
    geometries are (Multi)LineStrings, they are converted into
    (Multi)Polygons.

    Args:
        gdf: GeoDataFrame with geometries that are either (Multi)LineString
            or (Multi)Polygon.

    Returns:
        polygon_gdf: Transformation of gdf where all geometries are
            (Multi)Polygon.
    """

    polygon_gdf = gdf.copy()

    for index, row in polygon_gdf.iterrows():

        # Extract geometry
        geom = row['geometry']

        # Check if geometry is not Polygon
        if isinstance(geom, (LineString, MultiLineString)):

            # Create new geometry that is a Polygon
            new_geom = convert_linestring_to_polygon(geom)

            # Store new Polygon geometry in row
            # We use .loc[index, 'geometry'] to change value

            if isinstance(geom, MultiLineString):
                polygon_gdf.loc[index, 'geometry'] = new_geom.values
            else:
                polygon_gdf.loc[index, 'geometry'] = new_geom

    return polygon_gdf

def identify_polygon_from_map_no(uac_gdf, map_no):
    """Identifies which polygon is most likely associated with map_no

    Finds all polygons associated with a specific Map No and computes distances
    from their centroids to the centroid of the MultiPolygon merging all these
    polygons. The polygon with the shortest distance is selected.

    Args:
        uac_gdf: GeoDataFrame with polygon geometries. It is assumed
            to be the data for unauthorized colonies.
        map_no: Map_No for which there are duplicate entries in
            polygon_gdf.

    Returns:
        polygons_to_keep_delete: Python distionary with two keys: "keep"
            and "delete". "keep" has as its value the index of the
            polygon to keep while "delete" has as its value a list
            of indices for all rows to delete.
    """

    # Create a subset of uac_gdf only consisting of rows
    # where "Map_No" equals map_no
    uac_gdf_subset =  uac_gdf[uac_gdf['MAP_NO'] == map_no]

    # Take the union of all geometries in the subset
    merged_multipolygon_map_no = uac_gdf_subset.unary_union

    # Compute the centroid of the merged MultiPolygon
    map_no_centroid = merged_multipolygon_map_no.centroid

    # Initialize dictionary with centroid distances
    # Key is index and value is distance
    centroid_distances = dict()

    # Iterate across all rows in subset matching Map_No
    for idx, row in uac_gdf_subset.iterrows():
        # Add centroid distances to centroid_distance dictionary
        centroid_distances[row['index']] = row['geometry'].\
            centroid.distance(map_no_centroid)

    # Calculate the minimum distance of all centroid distances
    min_distance = min(centroid_distances.values())

    # Initialize dictionary with polygons to keep and delete
    polygons_to_keep_delete = dict()

    # Initialize delete key with a value of an empty list
    polygons_to_keep_delete['delete'] = []

    # If index corresponds to shortest distance designate as "keep"
    # else place the index list for the "delete" key
    for index in centroid_distances.keys():
        if centroid_distances[index] == min_distance:
            polygons_to_keep_delete['keep'] = index
        else:
            polygons_to_keep_delete['delete'].append(index)

    return polygons_to_keep_delete

## RESULTS for identify_polygon_from_map_no
### Checking random 5% of duplicate map numbers
#* 733: Correct!
#* 497: Correct!
#* 348: Correct!
#* 1436: Correct!
#* 1572: Correct!
#* 865: Correct!
#* 974: Correct!
#* 738: Correct!
#* 990: Correct!
#* 224: Not Correct!
