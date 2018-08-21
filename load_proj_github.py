'''
6.Trimesh : 
1.trimesh.load_mesh	(load.py in io folder)
2.centroid	(base.py)	?
3.bounds	(bounds.py)	
4.apply_translation	(base.py)

'''


#1.trimesh.load_mesh	(load.py in io folder)
@_log_time
def load_mesh(file_obj, file_type=None, **kwargs):
    """
    Load a mesh file into a Trimesh object

    Parameters
    ---------
    file_obj:  str or file-like object
    file_type: str representing file type (eg: 'stl')
    kwargs:    passed to Trimesh constructor

    Returns:
    ----------
    mesh: Trimesh object, or a list of Trimesh objects
          depending on the file format.

    """
    # turn a string into a file obj and type
    (file_obj,
     file_type,
     metadata,
     opened) = _parse_file_args(file_obj, file_type)

    # make sure we keep passed kwargs to loader
    # but also make sure loader keys override passed keys
    results = mesh_loaders[file_type](file_obj,
                                      file_type=file_type)

    if util.is_file(file_obj):
        file_obj.close()

    log.debug('loaded mesh using %s',
              mesh_loaders[file_type].__name__)

    if not isinstance(results, list):
        results = [results]

    loaded = []
    for result in results:
        kwargs.update(result)
        loaded.append(load_kwargs(kwargs))
        loaded[-1].metadata.update(metadata)
    if len(loaded) == 1:
        loaded = loaded[0]

    if opened:
        file_obj.close()

    return loaded
    
#--------------------------

def load_kwargs(*args, **kwargs):
    """
    Load geometry from a properly formatted dict or kwargs

    """
    def handle_scene():
        """
        Load a scene from our kwargs:

        class:      Scene
        geometry:   dict, name: Trimesh kwargs
        graph:      list of dict, kwargs for scene.graph.update
        base_frame: str, base frame of graph
        """
        scene = Scene()
        scene.geometry.update({k: load_kwargs(v) for
                               k, v in kwargs['geometry'].items()})
        for k in kwargs['graph']:
            if isinstance(k, dict):
                scene.graph.update(**k)
            elif util.is_sequence(k) and len(k) == 3:
                scene.graph.update(k[1], k[0], **k[2])
        if 'base_frame' in kwargs:
            scene.graph.base_frame = kwargs['base_frame']
        if 'metadata' in kwargs:
            scene.metadata.update(kwargs['metadata'])

        return scene

    def handle_trimesh_kwargs():
        if (isinstance(kwargs['vertices'], dict) or
                isinstance(kwargs['faces'], dict)):
            return Trimesh(**misc.load_dict(kwargs))
        else:
            return Trimesh(**kwargs)

    def handle_trimesh_export():
        data, file_type = kwargs['data'], kwargs['file_type']
        if not isinstance(data, dict):
            data = util.wrap_as_stream(data)
        k = mesh_loaders[file_type](data,
                                    file_type=file_type)
        return Trimesh(**k)

    # if we've been passed a single dict instead of kwargs
    # substitute the dict for kwargs
    if (len(kwargs) == 0 and
        len(args) == 1 and
            isinstance(args[0], dict)):
        kwargs = args[0]

    handlers = {handle_scene: ('graph', 'geometry'),
                handle_trimesh_kwargs: ('vertices', 'faces'),
                handle_trimesh_export: ('file_type', 'data')}

    handler = None
    for k, v in handlers.items():
        if all(i in kwargs for i in v):
            handler = k
            break
    if handler is None:
        raise ValueError('unable to determine type!')

    return handler()


def _parse_file_args(file_obj, file_type, **kwargs):
    """
    Given a file_obj and a file_type, try to turn them into a file-like object
    and a lowercase string of file type

    Parameters
    -----------
    file_obj:  str: if string represents a file path, returns
                    -------------------------------------------
                    file_obj:   an 'rb' opened file object of the path
                    file_type:  the extension from the file path

               str: if string is NOT a path, but has JSON-like special characters
                    -------------------------------------------
                    file_obj:   the same string passed as file_obj
                    file_type:  set to 'json'

               str: string is not an existing path or a JSON-like object
                    -------------------------------------------
                    ValueError will be raised as we can't do anything with input

               file like object: we cannot grab information on file_type automatically
                    -------------------------------------------
                    ValueError will be raised if file_type is None
                    file_obj:  same as input
                    file_type: same as input

               other object: like a shapely.geometry.Polygon, etc:
                    -------------------------------------------
                    file_obj:  same as input
                    file_type: if None initially, set to the class name
                               (in lower case), otherwise passed through

    file_type: str, type of file and handled according to above

    Returns
    -----------
    file_obj:  loadable object
    file_type: str, lower case of the type of file (eg 'stl', 'dae', etc)
    metadata:  dict, any metadata
    opened:    bool, did we open the file or not
    """
    metadata = {}
    opened = False
    if ('metadata' in kwargs and
            isinstance(kwargs['metadata'], dict)):
        metadata.update(kwargs['metadata'])

    if util.is_file(file_obj) and file_type is None:
        raise ValueError(
            'File type must be specified when passing file objects!')
    if util.is_string(file_obj):
        try:
            # os.path.isfile will return False incorrectly
            # if we don't give it an absolute path
            file_path = os.path.expanduser(file_obj)
            file_path = os.path.abspath(file_path)
            exists = os.path.isfile(file_path)
        except BaseException:
            exists = False

        if exists:
            metadata['file_path'] = file_path
            metadata['file_name'] = os.path.basename(file_obj)
            # if file_obj is a path that exists use extension as file_type
            if file_type is None:
                file_type = util.split_extension(file_path,
                                                 special=['tar.gz', 'tar.bz2'])
            file_obj = open(file_path, 'rb')
            opened = True
        else:
            if file_type is not None:
                return file_obj, file_type, metadata
            elif '{' in file_obj:
                # if a dict bracket is in the string, its probably a straight
                # JSON
                file_type = 'json'
            else:
                raise ValueError(
                    'File object passed as string that is not a file!')

    if file_type is None:
        file_type = file_obj.__class__.__name__

    if util.is_string(file_type) and '.' in file_type:
        # if someone has passed the whole filename as the file_type
        # use the file extension as the file_type
        if 'file_path' not in metadata:
            metadata['file_path'] = file_type
        metadata['file_name'] = os.path.basename(file_type)
        file_type = util.split_extension(file_type)
    file_type = file_type.lower()
    return file_obj, file_type, metadata, opened


compressed_loaders = {'zip': load_compressed,
                      'tar.bz2': load_compressed,
                      'tar.gz': load_compressed}
mesh_loaders = {}
# assimp has a lot of loaders, but they are all quite slow
# so we load them first and replace them with native loaders if possible
mesh_loaders.update(_assimp_loaders)
mesh_loaders.update(_stl_loaders)
mesh_loaders.update(_ctm_loaders)
mesh_loaders.update(_misc_loaders)
mesh_loaders.update(_ply_loaders)
mesh_loaders.update(_xml_loaders)
mesh_loaders.update(_obj_loaders)
mesh_loaders.update(_gltf_loaders)
mesh_loaders.update(_three_loaders)












    

#2.centroid	(base.py)
@util.cache_decorator
def centroid(self):
        """
        The point in space which is the average of the triangle centroids
        weighted by the area of each triangle.

        This will be valid even for non- watertight meshes,
        unlike self.center_mass

        Returns
        ----------
        centroid: (3,) float, the average vertex
        """

        # use the centroid of each triangle weighted by
        # the area of the triangle to find the overall centroid
        centroid = np.average(self.triangles_center,
                              axis=0,
                              weights=self.area_faces)
        centroid.flags.writeable = False
        return centroid
        
#3.bounds	(base.py)
@util.cache_decorator
def bounds(self):
        """
        The axis aligned bounds of the mesh.

        Returns
        -----------
        bounds: (2,3) float, bounding box with [min, max] coordinates
        """
        # we use triangles instead of faces because
        # if there is an unused vertex it will screw up bounds
        in_mesh = self.triangles.reshape((-1, 3))
        bounds = np.vstack((in_mesh.min(axis=0),
                            in_mesh.max(axis=0)))
        bounds.flags.writeable = False
        return bounds
        
#4.apply_translation	(base.py)
def apply_translation(self, translation):
        """
        Translate the current mesh.

        Parameters
        ----------
        translation: (3,) float, translation in XYZ
        """
        translation = np.asanyarray(translation, dtype=np.float64)
        if translation.shape != (3,):
            raise ValueError('Translation must be (3,)!')

        matrix = np.eye(4)
        matrix[:3, 3] = translation
        self.apply_transform(matrix)