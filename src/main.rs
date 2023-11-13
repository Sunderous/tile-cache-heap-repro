use bitflags::bitflags;
use fastlz_sys::{fastlz_compress, fastlz_decompress};
use recastnavigation_sys::*;
use std::{
    fs,
    os::raw::{c_uchar, c_void},
};

pub fn main() {}

#[repr(C)]
struct TileCacheSetHeader {
    magic: i32,
    version: i32,
    num_tiles: i32,
    mesh_params: dtNavMeshParams,
    cache_params: dtTileCacheParams,
}

#[repr(C)]
struct TileCacheTileHeader {
    tile_ref: dtCompressedTileRef,
    data_size: i32,
}

// The flags for points in a "straight path".
// Note: recastnavigation-sys generates either i32 or u32 for enums, but dtStraightPathFlags are
// generally taken as u8.
bitflags! {
    #[repr(transparent)]
    pub struct DtStraightPathFlags: u8 {
        const START = dtStraightPathFlags_DT_STRAIGHTPATH_START as _;
        const END = dtStraightPathFlags_DT_STRAIGHTPATH_END as _;
        const OFFMESH_CONNECTION = dtStraightPathFlags_DT_STRAIGHTPATH_OFFMESH_CONNECTION as _;
    }
}

// extern "C" fn max_compressed_size(_object_ptr: *mut std::ffi::c_void, buffer_size: i32) -> i32 {
//     buffer_size
// }

// extern "C" fn compress(
//     _object_ptr: *mut std::ffi::c_void,
//     buffer: *const u8,
//     buffer_size: i32,
//     compressed: *mut u8,
//     max_compressed_size: i32,
//     compressed_size: *mut i32,
// ) -> u32 {
//     assert!(
//         buffer_size <= max_compressed_size,
//         "\n\nleft: {}\nright: {}",
//         buffer_size,
//         max_compressed_size
//     );

//     let buffer_slice = unsafe { std::slice::from_raw_parts(buffer, buffer_size as usize) };

//     unsafe { *compressed_size = buffer_size };

//     let compressed_slice =
//         unsafe { std::slice::from_raw_parts_mut(compressed, *compressed_size as usize) };

//     compressed_slice.copy_from_slice(buffer_slice);

//     DT_SUCCESS
// }

// extern "C" fn decompress(
//     object_ptr: *mut std::ffi::c_void,
//     compressed: *const u8,
//     compressed_size: i32,
//     buffer: *mut u8,
//     max_buffer_size: i32,
//     buffer_size: *mut i32,
// ) -> u32 {
//     // Since compress just copies the source to destination, decompress is
//     // the exact same.
//     compress(
//         object_ptr,
//         compressed,
//         compressed_size,
//         buffer,
//         max_buffer_size,
//         buffer_size,
//     )
// }

extern "C" fn max_compressed_size(_object_ptr: *mut std::ffi::c_void, buffer_size: i32) -> i32 {
    return (buffer_size as f32 * 1.05f32) as i32;
}

extern "C" fn compress(
    _object_ptr: *mut std::ffi::c_void,
    buffer: *const u8,
    buffer_size: i32,
    compressed: *mut u8,
    _max_compressed_size: i32,
    compressed_size: *mut i32,
) -> u32 {
    unsafe {
        *compressed_size = fastlz_compress(
            buffer as *const c_void,
            buffer_size,
            compressed as *mut c_void,
        );
    }

    return DT_SUCCESS;
}

extern "C" fn decompress(
    _object_ptr: *mut std::ffi::c_void,
    compressed: *const u8,
    compressed_size: i32,
    buffer: *mut u8,
    max_buffer_size: i32,
    buffer_size: *mut i32,
) -> u32 {
    unsafe {
        *buffer_size = fastlz_decompress(
            compressed as *const c_void,
            compressed_size,
            buffer as *mut c_void,
            max_buffer_size,
        );

        if *buffer_size < 0 {
            return DT_FAILURE;
        }

        return DT_SUCCESS;
    }
}

#[test]
fn detour_tile_cache_simple_caching() {
    let cache_params = dtTileCacheParams {
        orig: [0.0, 0.0, 0.0],
        cs: 1.0,
        ch: 1.0,
        width: 5,
        height: 5,
        walkableHeight: 1.0,
        walkableRadius: 1.0,
        walkableClimb: 1.0,
        maxSimplificationError: 0.01,
        maxTiles: 1000,
        maxObstacles: 10,
    };

    let alloc = unsafe { CreateDefaultTileCacheAlloc() };

    let forwarded_compressor = unsafe {
        CreateForwardedTileCacheCompressor(
            std::ptr::null_mut(),
            Some(max_compressed_size),
            Some(compress),
            Some(decompress),
        )
    };

    extern "C" fn set_poly_flags(
        _: *mut std::ffi::c_void,
        params: *mut dtNavMeshCreateParams,
        _areas: *mut u8,
        flags: *mut u16,
    ) {
        let params = unsafe { &*params };
        let flags = unsafe { std::slice::from_raw_parts_mut(flags, params.polyCount as usize) };

        flags.fill(1);
    }
    let forwarded_mesh_process =
        unsafe { CreateForwardedTileCacheMeshProcess(std::ptr::null_mut(), Some(set_poly_flags)) };

    let tile_cache = unsafe { &mut *dtAllocTileCache() };
    unsafe {
        tile_cache.init(
            &cache_params,
            alloc,
            forwarded_compressor,
            forwarded_mesh_process,
        )
    };

    let nav_mesh = unsafe { &mut *dtAllocNavMesh() };
    let nav_mesh_params = dtNavMeshParams {
        orig: [0.0, 0.0, 0.0],
        tileWidth: 5.0,
        tileHeight: 5.0,
        maxTiles: 1000,
        maxPolys: 10,
    };
    unsafe { nav_mesh.init(&nav_mesh_params) };

    for _ in 0..10 {
        let mut up_to_date = false;
        unsafe { tile_cache.update(1.0, nav_mesh, &mut up_to_date) };
        assert!(up_to_date);
    }

    let mut header = dtTileCacheLayerHeader {
        magic: DT_TILECACHE_MAGIC,
        version: DT_TILECACHE_VERSION,
        tx: 0,
        ty: 0,
        tlayer: 0,
        bmin: [0.0, 1.0, 0.0],
        bmax: [5.0, 1.0, 5.0],
        width: 5 as u8,
        height: 5 as u8,
        minx: 0,
        maxx: 4,
        miny: 0,
        maxy: 4,
        hmin: 1,
        hmax: 1,
    };

    const N: u8 = 255;

    let heights = [
        N, N, 0, N, N, //
        N, N, 0, N, N, //
        N, N, 0, N, N, //
        N, N, 0, N, N, //
        0, 0, 0, 0, 0, //
    ];

    const W: u8 = DT_TILECACHE_WALKABLE_AREA;

    let areas = [
        0, 0, W, 0, 0, //
        0, 0, W, 0, 0, //
        0, 0, W, 0, 0, //
        0, 0, W, 0, 0, //
        W, W, W, W, W, //
    ];

    // Neighbour connectivity.
    let cons = [
        0, 0, 2, 0, 0, //
        0, 0, 10, 0, 0, //
        0, 0, 10, 0, 0, //
        0, 0, 10, 0, 0, //
        4, 5, 13, 5, 1, //
    ];

    let mut data: *mut u8 = std::ptr::null_mut();
    let mut data_size: i32 = 0;

    assert_eq!(
        unsafe {
            dtBuildTileCacheLayer(
                forwarded_compressor,
                &mut header,
                heights.as_ptr(),
                areas.as_ptr(),
                cons.as_ptr(),
                &mut data,
                &mut data_size,
            )
        },
        DT_SUCCESS
    );

    assert_eq!(
        unsafe {
            tile_cache.addTile(
                data,
                data_size,
                dtTileFlags_DT_TILE_FREE_DATA as u8,
                std::ptr::null_mut(),
            )
        },
        DT_SUCCESS
    );

    assert_eq!(
        unsafe { tile_cache.buildNavMeshTilesAt(0, 0, nav_mesh) },
        DT_SUCCESS
    );

    let query = unsafe { &mut *dtAllocNavMeshQuery() };
    assert_eq!(unsafe { query.init(nav_mesh, 10) }, DT_SUCCESS);

    let query_filter = dtQueryFilter {
        m_areaCost: [1.0; 64],
        m_includeFlags: 0xffff,
        m_excludeFlags: 0,
    };

    let mut path = [0; 10];
    let mut path_count = 0;

    let start_point = [2.1, 1.0, 0.1];
    let end_point = [4.9, 1.0, 4.9];

    let mut start_point_ref = 0;
    assert_eq!(
        unsafe {
            query.findNearestPoly(
                start_point.as_ptr(),
                [0.1, 100.0, 0.1].as_ptr(),
                &query_filter,
                &mut start_point_ref,
                std::ptr::null_mut(),
            )
        },
        DT_SUCCESS
    );
    assert_ne!(start_point_ref, 0);

    let mut end_point_ref = 0;
    assert_eq!(
        unsafe {
            query.findNearestPoly(
                end_point.as_ptr(),
                [0.1, 100.0, 0.1].as_ptr(),
                &query_filter,
                &mut end_point_ref,
                std::ptr::null_mut(),
            )
        },
        DT_SUCCESS
    );
    assert_ne!(end_point_ref, 0);

    assert_eq!(
        unsafe {
            query.findPath(
                start_point_ref,
                end_point_ref,
                start_point.as_ptr(),
                end_point.as_ptr(),
                &query_filter,
                path.as_mut_ptr(),
                &mut path_count,
                path.len() as i32,
            )
        },
        DT_SUCCESS
    );

    assert_eq!(
        path.iter()
            .map(|polyref| if *polyref == 0 {
                None
            } else {
                Some(polyref & 0b11111111111111)
            })
            .collect::<Vec<_>>(),
        [
            Some(1),
            Some(3),
            Some(0),
            None,
            None,
            None,
            None,
            None,
            None,
            None
        ]
    );

    unsafe { dtFreeNavMeshQuery(query) };
    unsafe { dtFreeNavMesh(nav_mesh) };
    unsafe { dtFreeTileCache(tile_cache) };
    unsafe { DeleteTileCacheMeshProcess(forwarded_mesh_process) };
    unsafe { DeleteTileCacheCompressor(forwarded_compressor) };
    unsafe { DeleteTileCacheAlloc(alloc) };
}

#[test]
fn test_loading_tile_cache() {
    unsafe {
        const SAMPLE_POLYFLAGS_DISABLED: u16 = 0x10; // Disabled polygon
        const SAMPLE_POLYFLAGS_ALL: u16 = 0xffff; // All abilities

        let (mesh, cache, alloc, compressor, process) =
            load_nav_tile_cache("./nav_meshes/sample_tile_cache.dtcache");

        let nav_query = &mut *dtAllocNavMeshQuery();
        nav_query.init(mesh, 2048);

        let mut filter = dtQueryFilter::new();
        filter.m_includeFlags = SAMPLE_POLYFLAGS_ALL ^ SAMPLE_POLYFLAGS_DISABLED;
        filter.m_excludeFlags = 0;

        let start_pos = [-18.777670, 1.495857, -15.935596];
        let end_pos = [-56.686096, 1.257904, 53.642925];

        let start_ref = find_poly_ref(&nav_query, &start_pos);
        let end_ref = find_poly_ref(&nav_query, &end_pos);

        println!("Start Ref: {}", start_ref);
        println!("End Ref: {}", end_ref);

        let mut poly_path = Vec::<u32>::with_capacity(256);
        let mut poly_path_count: i32 = 0;

        nav_query.findPath(
            start_ref,
            end_ref,
            start_pos.as_ptr(),
            end_pos.as_ptr(),
            &filter,
            poly_path.as_mut_ptr() as *mut dtPolyRef,
            (&mut poly_path_count) as *mut i32,
            256,
        );

        poly_path.set_len(poly_path_count.try_into().unwrap());

        let mut straight_path_points = Vec::<[f32; 3]>::with_capacity(256);
        let mut straight_path_flags = [0; 256];
        let mut straight_path_polys = [0; 256];
        let mut straight_path_count = 0;

        let result = nav_query.findStraightPath(
            start_pos.as_ptr(),
            end_pos.as_ptr(),
            poly_path.as_ptr(),
            poly_path.len().try_into().unwrap(),
            straight_path_points.as_mut_ptr() as *mut f32,
            straight_path_flags.as_mut_ptr() as _,
            straight_path_polys.as_mut_ptr() as *mut dtPolyRef,
            &mut straight_path_count,
            256,
            0,
        );

        straight_path_points.set_len(straight_path_count.try_into().unwrap());

        println!("Result Value: {}", result);
        println!("Result Path Count: {}", straight_path_count);
        println!("Path: {:?}", straight_path_points);

        assert!(straight_path_points.len() == 12);

        dtFreeNavMeshQuery(nav_query);
        dtFreeNavMesh(mesh);

        // TODO: Why does this corrupt the heap?
        dtFreeTileCache(cache);

        DeleteTileCacheMeshProcess(process);
        DeleteTileCacheCompressor(compressor);
        DeleteTileCacheAlloc(alloc);
    }
}

fn find_poly_ref(query: &dtNavMeshQuery, pos: &[f32; 3]) -> dtPolyRef {
    let extents = [0.1, 100.0, 0.1];

    let mut poly_ref: dtPolyRef = 0;

    let query_filter = dtQueryFilter {
        m_areaCost: [1.0; 64],
        m_includeFlags: 0xffff,
        m_excludeFlags: 0,
    };

    assert_eq!(
        unsafe {
            query.findNearestPoly(
                pos.as_ptr(),
                extents.as_ptr(),
                &query_filter,
                &mut poly_ref,
                std::ptr::null_mut(),
            )
        },
        DT_SUCCESS
    );

    poly_ref
}

#[repr(u8)]
enum SamplePolyAreas {
    SAMPLE_POLYAREA_GROUND = 0,
    SAMPLE_POLYAREA_WATER,
    SAMPLE_POLYAREA_ROAD,
    SAMPLE_POLYAREA_DOOR,
    SAMPLE_POLYAREA_GRASS,
    SAMPLE_POLYAREA_JUMP,
}

#[repr(u16)]
enum SamplePolyFlags {
    SAMPLE_POLYFLAGS_WALK = 0x01, // Ability to walk (ground, grass, road)
    SAMPLE_POLYFLAGS_SWIM = 0x02, // Ability to swim (water).
    SAMPLE_POLYFLAGS_DOOR = 0x04, // Ability to move through doors.
    SAMPLE_POLYFLAGS_JUMP = 0x08, // Ability to jump.
    SAMPLE_POLYFLAGS_DISABLED = 0x10, // Disabled polygon
    SAMPLE_POLYFLAGS_ALL = 0xffff, // All abilities.
}

extern "C" fn set_poly_flags(
    _: *mut std::ffi::c_void,
    params: *mut dtNavMeshCreateParams,
    areas: *mut u8,
    flags: *mut u16,
) {
    // Update poly flags from areas.
    unsafe {
        let mut i: usize = 0;
        let poly_count = (*params).polyCount;
        let poly_areas = std::slice::from_raw_parts_mut(areas, poly_count as usize);
        let poly_flags = std::slice::from_raw_parts_mut(flags, poly_count as usize);

        while i < poly_count as usize {
            if poly_areas[i] == DT_TILECACHE_WALKABLE_AREA {
                poly_areas[i] = SamplePolyAreas::SAMPLE_POLYAREA_GROUND as u8;
            }

            if poly_areas[i] == SamplePolyAreas::SAMPLE_POLYAREA_GROUND as u8
                || poly_areas[i] == SamplePolyAreas::SAMPLE_POLYAREA_GRASS as u8
                || poly_areas[i] == SamplePolyAreas::SAMPLE_POLYAREA_ROAD as u8
            {
                poly_flags[i] = SamplePolyFlags::SAMPLE_POLYFLAGS_WALK as u16;
            } else if poly_areas[i] == SamplePolyAreas::SAMPLE_POLYAREA_WATER as u8 {
                poly_flags[i] = SamplePolyFlags::SAMPLE_POLYFLAGS_SWIM as u16;
            } else if poly_areas[i] == SamplePolyAreas::SAMPLE_POLYAREA_DOOR as u8 {
                poly_flags[i] = SamplePolyFlags::SAMPLE_POLYFLAGS_WALK as u16
                    | SamplePolyFlags::SAMPLE_POLYFLAGS_DOOR as u16;
            }

            i += 1;
        }
    }
}

unsafe fn load_nav_tile_cache(
    path: &str,
) -> (
    &mut dtNavMesh,
    &mut dtTileCache,
    *mut dtTileCacheAlloc,
    *mut dtTileCacheCompressor,
    *mut dtTileCacheMeshProcess,
) {
    const TILECACHESET_MAGIC: i32 =
        ('T' as i32) << 24 | ('S' as i32) << 16 | ('E' as i32) << 8 | ('T' as i32); //'TSET';
    const TILECACHESET_VERSION: i32 = 1;

    let mut data = fs::read(path).expect("Unable to read nav mesh file");
    let header_bytes: [u8; 92] = data[..92].try_into().unwrap();
    let header = std::mem::transmute::<[u8; 92], TileCacheSetHeader>(header_bytes);
    println!("Header - magic: {}", header.magic);
    println!("Header - version: {}", header.version);
    println!("Header - numtiles: {}", header.num_tiles);
    println!("Params - tileWidth: {}", header.mesh_params.tileWidth);
    println!("Params - tileHeight: {}", header.mesh_params.tileHeight);
    println!("Params - maxTiles: {}", header.mesh_params.maxTiles);
    println!("Params - maxPolys: {}", header.mesh_params.maxPolys);

    if header.magic != TILECACHESET_MAGIC {
        println!("FAILED MAGIC CHECK");
    }
    if header.version != TILECACHESET_VERSION {
        println!("FAILED VERSION CHECK");
    }

    let mesh = &mut *dtAllocNavMesh();

    let mut status = mesh.init(&header.mesh_params);

    // TODO: How to get access/re-implement the dtStatusFailed() function?
    if (status & DT_FAILURE) != 0 {
        println!("FAILED TO INITIALIZE MESH WITH PARAMS");
    }

    let cache = &mut *dtAllocTileCache();

    let default_alloc = unsafe { CreateDefaultTileCacheAlloc() };

    let forwarded_compressor = unsafe {
        CreateForwardedTileCacheCompressor(
            std::ptr::null_mut(),
            Some(max_compressed_size),
            Some(compress),
            Some(decompress),
        )
    };

    let forwarded_mesh_process =
        unsafe { CreateForwardedTileCacheMeshProcess(std::ptr::null_mut(), Some(set_poly_flags)) };

    status = cache.init(
        &header.cache_params,
        default_alloc,
        forwarded_compressor,
        forwarded_mesh_process,
    );

    if (status & DT_FAILURE) != 0 {
        println!("FAILED TO INITIALIZE CACHE WITH PARAMS");
    }

    let data_bytes = &mut data.as_mut_slice()[92..];
    let mut position = 0;

    for _ in 0..header.num_tiles {
        let tile_header_bytes: [u8; 8] = data_bytes[position..position + 8].try_into().unwrap();
        let tile_header = std::mem::transmute::<[u8; 8], TileCacheTileHeader>(tile_header_bytes);

        position += 8;

        let data = &mut data_bytes[position..position + (tile_header.data_size as usize)];

        let mut tile: dtCompressedTileRef = 0;

        cache.addTile(
            data.as_mut_ptr(),
            tile_header.data_size,
            dtCompressedTileFlags_DT_COMPRESSEDTILE_FREE_DATA as c_uchar,
            &mut tile as *mut dtCompressedTileRef,
        );

        if tile > 0 {
            cache.buildNavMeshTile(tile, mesh);
        }

        position += tile_header.data_size as usize;
    }

    return (
        mesh,
        cache,
        default_alloc,
        forwarded_compressor,
        forwarded_mesh_process,
    );
}
