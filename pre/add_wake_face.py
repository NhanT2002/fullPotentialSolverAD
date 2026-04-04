#!/usr/bin/env python3

import argparse
import shutil
from pathlib import Path

import h5py
import numpy as np


def set_cgns_name(node, name):
	name_bytes = name.encode("ascii")[:32].ljust(32, b"\x00")
	value = np.frombuffer(name_bytes, dtype="S1")
	if "name" in node.attrs:
		node.attrs["name"] = value
	else:
		node.attrs.create("name", value)


def get_cgns_name(node):
	if "name" not in node.attrs:
		return ""
	raw = node.attrs["name"]
	if isinstance(raw, bytes):
		return raw.decode("ascii", errors="ignore").rstrip("\x00")
	if isinstance(raw, np.ndarray):
		try:
			return b"".join(raw.tolist()).decode("ascii", errors="ignore").rstrip("\x00")
		except TypeError:
			return "".join(str(v, "ascii") for v in raw).rstrip("\x00")
	return str(raw)


def replace_data_dataset(group, data):
	if " data" in group:
		del group[" data"]
	group.create_dataset(" data", data=data)


def set_cgns_core_attrs(node, label, type_code):
	node.attrs["flags"] = np.array([1], dtype=np.int32)
	node.attrs["label"] = np.bytes_(label)
	node.attrs["type"] = np.bytes_(type_code)


def ensure_cgns_group(parent, name, label, type_code):
	if name in parent:
		group = parent[name]
	else:
		group = parent.create_group(name)
	set_cgns_core_attrs(group, label, type_code)
	set_cgns_name(group, name)
	return group


def ascii_i8(text):
	return np.frombuffer(text.encode("ascii"), dtype=np.int8)


def ensure_bc_child_group(parent, name):
	if name in parent:
		return parent[name]
	return parent.create_group(name)


def ensure_zonebc_entry(zone_bc, name, template=None):
	if name in zone_bc:
		bc = zone_bc[name]
	else:
		if template is not None:
			zone_bc.copy(template, name)
			bc = zone_bc[name]
		else:
			bc = zone_bc.create_group(name)
			replace_data_dataset(bc, ascii_i8("FamilySpecified"))

	set_cgns_name(bc, name)

	point_range = ensure_bc_child_group(bc, "PointRange")
	grid_location = ensure_bc_child_group(bc, "GridLocation")
	family_name = ensure_bc_child_group(bc, "FamilyName")

	set_cgns_name(point_range, "PointRange")
	set_cgns_name(grid_location, "GridLocation")
	set_cgns_name(family_name, "FamilyName")

	if " data" not in grid_location:
		replace_data_dataset(grid_location, ascii_i8("EdgeCenter"))

	return bc


def ensure_base_wake_families(base, family_names, enabled=True):
	if not enabled:
		for name in ["wake", "wake_lower", "wake_upper"]:
			if name in base:
				del base[name]
		return

	template_name = "farfield" if "farfield" in base else ("wall" if "wall" in base else None)

	for fam_name in family_names:
		if fam_name in base:
			wake_family = base[fam_name]
		else:
			if template_name is not None:
				base.copy(template_name, fam_name)
				wake_family = base[fam_name]
			else:
				wake_family = base.create_group(fam_name)

		set_cgns_name(wake_family, fam_name)

		if template_name is not None:
			template_family = base[template_name]
			for child_name in template_family.keys():
				if child_name not in wake_family:
					wake_family.copy(template_family[child_name], child_name)

		for child_name in ["FamBC", "Fam_Descr_Name", "FamBC_TypeId", "FamBC_TypeName", "FamBC_UserId", "FamBC_UserName"]:
			if child_name not in wake_family:
				wake_family.create_group(child_name)
			set_cgns_name(wake_family[child_name], child_name)

		# Keep wake family names, but map BC type to wall for Pointwise compatibility.
		if "wall" in base:
			wall_family = base["wall"]
			for key in ["FamBC", "FamBC_TypeId", "FamBC_TypeName", "FamBC_UserId"]:
				if key in wall_family and " data" in wall_family[key]:
					replace_data_dataset(wake_family[key], wall_family[key][" data"][:])
				elif key == "FamBC":
					replace_data_dataset(wake_family[key], ascii_i8("BCWallInviscid"))
				elif key == "FamBC_TypeId":
					replace_data_dataset(wake_family[key], ascii_i8("21"))
				elif key == "FamBC_TypeName":
					replace_data_dataset(wake_family[key], ascii_i8("Wall Inviscid"))
				else:
					replace_data_dataset(wake_family[key], ascii_i8("1"))
		else:
			replace_data_dataset(wake_family["FamBC"], ascii_i8("BCWallInviscid"))
			replace_data_dataset(wake_family["FamBC_TypeId"], ascii_i8("21"))
			replace_data_dataset(wake_family["FamBC_TypeName"], ascii_i8("Wall Inviscid"))
			replace_data_dataset(wake_family["FamBC_UserId"], ascii_i8("1"))

		replace_data_dataset(wake_family["Fam_Descr_Name"], ascii_i8(fam_name))
		replace_data_dataset(wake_family["FamBC_UserName"], ascii_i8(fam_name))

	if "wake" in base and "wake" not in family_names:
		del base["wake"]


def update_zonebc_tags(zone, farfield_range, wall_range, wake_ranges=None):
	zone_bc = zone["ZoneBC"] if "ZoneBC" in zone else zone.create_group("ZoneBC")
	set_cgns_name(zone_bc, "ZoneBC")

	template = None
	for key in ("farfield", "wall", "wake_lower", "wake_upper", "wake"):
		if key in zone_bc:
			template = zone_bc[key]
			break

	entries = {
		"farfield": farfield_range,
		"wall": wall_range,
	}
	if wake_ranges is not None:
		entries.update(wake_ranges)

	for bc_name, elem_range in entries.items():
		bc = ensure_zonebc_entry(zone_bc, bc_name, template=template)
		replace_data_dataset(bc["PointRange"], np.array(elem_range, dtype=np.int64).reshape(2, 1))
		replace_data_dataset(bc["GridLocation"], ascii_i8("EdgeCenter"))
		replace_data_dataset(bc["FamilyName"], ascii_i8(bc_name))

	if wake_ranges is None:
		for name in ["wake", "wake_lower", "wake_upper"]:
			if name in zone_bc:
				del zone_bc[name]
	else:
		for name in ["wake", "wake_lower", "wake_upper"]:
			if name not in entries and name in zone_bc:
				del zone_bc[name]


def find_base_and_zone(handle):
	excluded = {"CGNSLibraryVersion", " format", " hdf5version"}
	base_name = None
	for name in handle.keys():
		if name not in excluded and isinstance(handle[name], h5py.Group):
			base_name = name
			break
	if base_name is None:
		raise RuntimeError("Could not find CGNS base node")

	base = handle[base_name]
	zone_name = None
	for name in base.keys():
		if isinstance(base[name], h5py.Group) and "GridCoordinates" in base[name]:
			zone_name = name
			break
	if zone_name is None:
		raise RuntimeError("Could not find CGNS zone node")

	return base, zone_name, base[zone_name]


def find_volume_section(zone):
	candidates = []
	for key, grp in zone.items():
		if not isinstance(grp, h5py.Group):
			continue
		if "ElementConnectivity" not in grp or " data" not in grp["ElementConnectivity"]:
			continue
		if " data" not in grp:
			continue
		sec_data = np.array(grp[" data"][:]).reshape(-1)
		if sec_data.size == 0:
			continue
		elem_type = int(sec_data[0])
		if elem_type in (5, 7):
			n_conn = int(grp["ElementConnectivity"][" data"].size)
			candidates.append((n_conn, key, elem_type))

	if not candidates:
		raise RuntimeError("Could not find TRI_3 or QUAD_4 volume section")

	_, sec_name, elem_type = max(candidates)
	nodes_per_elem = 3 if elem_type == 5 else 4
	return sec_name, nodes_per_elem


def find_named_bar_section(zone, target_name):
	target = target_name.lower()
	for key, grp in zone.items():
		if not isinstance(grp, h5py.Group):
			continue
		if "ElementConnectivity" not in grp or " data" not in grp["ElementConnectivity"]:
			continue
		if " data" not in grp:
			continue
		sec_data = np.array(grp[" data"][:]).reshape(-1)
		if sec_data.size == 0 or int(sec_data[0]) != 3:
			continue
		node_name = get_cgns_name(grp).lower()
		if key.lower() == target or node_name == target:
			return key
	return None


def read_airfoil_te(airfoil_file):
	coords = []
	with open(airfoil_file, "r", encoding="utf-8") as handle:
		for line in handle:
			stripped = line.strip()
			if not stripped:
				continue
			parts = stripped.replace(",", " ").split()
			if len(parts) < 2:
				continue
			try:
				x = float(parts[0])
				y = float(parts[1])
			except ValueError:
				continue
			coords.append((x, y))

	if len(coords) < 2:
		raise RuntimeError(f"No usable (x,y) points found in airfoil file: {airfoil_file}")

	arr = np.array(coords)
	x = arr[:, 0]
	y = arr[:, 1]
	x_te = float(np.max(x))

	x_span = max(1.0, float(np.max(x) - np.min(x)))
	tol = 1e-8 * x_span
	mask = np.abs(x - x_te) <= tol
	y_te = float(np.mean(y[mask])) if np.any(mask) else 0.0
	return x_te, y_te


def detect_te_from_wall(x, y, wall_edges):
	wall_nodes = np.unique(wall_edges.reshape(-1))
	x_wall = x[wall_nodes - 1]
	y_wall = y[wall_nodes - 1]

	x_te = float(np.max(x_wall))
	x_span = max(1.0, float(np.max(x_wall) - np.min(x_wall)))
	mask = np.abs(x_wall - x_te) <= 1e-8 * x_span
	y_te = float(np.mean(y_wall[mask])) if np.any(mask) else 0.0
	return x_te, y_te


def build_edge_to_elems(elem_conn):
	edge_to_elems = {}
	for elem_id, nodes in enumerate(elem_conn, start=1):
		for i, n1 in enumerate(nodes):
			n2 = nodes[(i + 1) % len(nodes)]
			key = (int(min(n1, n2)), int(max(n1, n2)))
			edge_to_elems.setdefault(key, []).append(elem_id)
	return edge_to_elems


def classify_wake_edges(edge_to_elems, x, y, x_te, y_te, y_tol, x_tol):
	wake_edges = []
	for (n1, n2), elems in edge_to_elems.items():
		y1 = y[n1 - 1]
		y2 = y[n2 - 1]
		if abs(y1 - y_te) > y_tol or abs(y2 - y_te) > y_tol:
			continue
		x_mid = 0.5 * (x[n1 - 1] + x[n2 - 1])
		if x_mid < x_te - x_tol:
			continue
		if len(elems) != 2:
			continue
		wake_edges.append((n1, n2, elems, x_mid))

	wake_edges.sort(key=lambda item: item[3])
	return wake_edges


def rewrite_boundary_with_duplicate_nodes(boundary_edges, y, y_te, y_tol, dup_node_map):
	updated = boundary_edges.copy()
	for i in range(updated.shape[0]):
		n1, n2 = int(updated[i, 0]), int(updated[i, 1])
		y_mid = 0.5 * (y[n1 - 1] + y[n2 - 1])
		if y_mid < y_te - y_tol:
			updated[i, 0] = dup_node_map.get(n1, n1)
			updated[i, 1] = dup_node_map.get(n2, n2)
	return updated


def update_section_connectivity(zone, section_name, conn_flat, elem_start, conn_dtype=None):
	section = zone[section_name]
	stored_conn_dtype = section["ElementConnectivity"][" data"].dtype
	stored_range_dtype = section["ElementRange"][" data"].dtype
	if conn_dtype is None:
		conn_dtype = stored_conn_dtype
	n_elem = conn_flat.size // 2
	elem_end = elem_start + n_elem - 1
	replace_data_dataset(section["ElementConnectivity"], conn_flat.astype(conn_dtype, copy=False))
	section["ElementRange"][" data"][...] = np.array([elem_start, elem_end], dtype=stored_range_dtype)
	return elem_end + 1


def write_wake_links_to_cgns(zone, wake_lower, wake_upper, wake_ranges, x_coords, y_coords, int_dtype):
	if "WakeLinkage" in zone:
		del zone["WakeLinkage"]

	udd = ensure_cgns_group(zone, "WakeLinkage", "UserDefinedData_t", "MT")

	n_pairs = len(wake_lower)
	pair_index = np.arange(1, n_pairs + 1, dtype=int_dtype)
	lower_cell = np.array([item[2] for item in wake_lower], dtype=int_dtype)
	upper_cell = np.array([item[2] for item in wake_upper], dtype=int_dtype)
	lower_n1 = np.array([item[0] for item in wake_lower], dtype=int_dtype)
	lower_n2 = np.array([item[1] for item in wake_lower], dtype=int_dtype)
	upper_n1 = np.array([item[0] for item in wake_upper], dtype=int_dtype)
	upper_n2 = np.array([item[1] for item in wake_upper], dtype=int_dtype)
	x_mid = np.array(
		[0.5 * (float(x_coords[item[0] - 1]) + float(x_coords[item[1] - 1])) for item in wake_upper],
		dtype=np.float64,
	)
	y_mid = np.array(
		[0.5 * (float(y_coords[item[0] - 1]) + float(y_coords[item[1] - 1])) for item in wake_upper],
		dtype=np.float64,
	)

	lower_face_elem_id = np.arange(
		int(wake_ranges["wake_lower"][0]), int(wake_ranges["wake_lower"][0]) + n_pairs, dtype=int_dtype
	)
	upper_face_elem_id = np.arange(
		int(wake_ranges["wake_upper"][0]), int(wake_ranges["wake_upper"][0]) + n_pairs, dtype=int_dtype
	)

	arrays = {
		"pairIndex": pair_index,
		"lowerCell": lower_cell,
		"upperCell": upper_cell,
		"lowerFaceElemId": lower_face_elem_id,
		"upperFaceElemId": upper_face_elem_id,
		"lowerNode1": lower_n1,
		"lowerNode2": lower_n2,
		"upperNode1": upper_n1,
		"upperNode2": upper_n2,
		"xMid": x_mid,
		"yMid": y_mid,
	}

	for name, values in arrays.items():
		type_code = "R8" if np.issubdtype(values.dtype, np.floating) else "I8"
		node = ensure_cgns_group(udd, name, "DataArray_t", type_code)
		replace_data_dataset(node, values)

	# Keep a one-line descriptor-like string for quick identification.
	desc = ensure_cgns_group(udd, "Description", "Descriptor_t", "C1")
	replace_data_dataset(desc, ascii_i8("Wake lower/upper linkage pairs"))


def split_wake(
	mesh_in,
	mesh_out,
	airfoil_file=None,
	y_tol=1e-10,
	x_tol=1e-12,
	overwrite=False,
	write_wake_boundary=True,
):
	src = Path(mesh_in)
	dst = Path(mesh_out)

	if not src.exists():
		raise FileNotFoundError(f"Input mesh not found: {src}")
	if dst.exists() and not overwrite:
		raise FileExistsError(f"Output mesh exists: {dst}. Use --overwrite to replace it.")

	if src.resolve() != dst.resolve():
		shutil.copy2(src, dst)

	with h5py.File(dst, "r+") as handle:
		base, _, zone = find_base_and_zone(handle)

		coords = zone["GridCoordinates"]
		x = coords["CoordinateX"][" data"][:]
		y = coords["CoordinateY"][" data"][:]
		z = coords["CoordinateZ"][" data"][:]

		wall_name = find_named_bar_section(zone, "wall")
		farfield_name = find_named_bar_section(zone, "farfield")
		if wall_name is None or farfield_name is None:
			raise RuntimeError("Could not find both wall and farfield BAR sections")

		vol_name, n_per_elem = find_volume_section(zone)
		vol_conn_dset = zone[vol_name]["ElementConnectivity"][" data"]
		vol_conn_dtype = vol_conn_dset.dtype
		vol_range_dtype = zone[vol_name]["ElementRange"][" data"].dtype
		vol_conn = vol_conn_dset[:]
		if vol_conn.size % n_per_elem != 0:
			raise RuntimeError("Volume connectivity size is inconsistent with element type")
		elem_conn = vol_conn.reshape(-1, n_per_elem)

		wall_conn_dtype = zone[wall_name]["ElementConnectivity"][" data"].dtype
		farfield_conn_dtype = zone[farfield_name]["ElementConnectivity"][" data"].dtype
		wall_edges = zone[wall_name]["ElementConnectivity"][" data"][:].reshape(-1, 2)
		farfield_edges = zone[farfield_name]["ElementConnectivity"][" data"][:].reshape(-1, 2)

		if airfoil_file:
			x_te, y_te = read_airfoil_te(airfoil_file)
		else:
			x_te, y_te = detect_te_from_wall(x, y, wall_edges)

		edge_to_elems = build_edge_to_elems(elem_conn)
		wake_edges = classify_wake_edges(edge_to_elems, x, y, x_te, y_te, y_tol, x_tol)

		if not wake_edges:
			raise RuntimeError(
				"No interior wake edges detected. Try increasing --y-tol or verify trailing-edge location."
			)

		elem_centroid_y = np.mean(y[elem_conn - 1], axis=1)

		wake_nodes = sorted(
			{n for (n1, n2, _, _) in wake_edges for n in (n1, n2)},
			key=lambda n: (x[n - 1], n),
		)
		dup_node_map = {}
		x_new = x.copy()
		y_new = y.copy()
		z_new = z.copy()
		next_node = x.size + 1
		for node in wake_nodes:
			dup_node_map[node] = next_node
			x_new = np.append(x_new, x[node - 1])
			y_new = np.append(y_new, y[node - 1])
			z_new = np.append(z_new, z[node - 1])
			next_node += 1

		lower_elems = set()
		wake_upper = []
		wake_lower = []

		for n1, n2, elems, _ in wake_edges:
			e1, e2 = elems[0], elems[1]
			if elem_centroid_y[e1 - 1] >= elem_centroid_y[e2 - 1]:
				upper_elem, lower_elem = e1, e2
			else:
				upper_elem, lower_elem = e2, e1

			lower_elems.add(lower_elem)

			if x[n1 - 1] <= x[n2 - 1]:
				up_n1, up_n2 = n1, n2
			else:
				up_n1, up_n2 = n2, n1

			wake_upper.append((up_n1, up_n2, upper_elem))
			wake_lower.append((dup_node_map[up_n1], dup_node_map[up_n2], lower_elem))

		elem_conn_new = elem_conn.copy()
		for elem in lower_elems:
			row = elem_conn_new[elem - 1]
			for i in range(row.size):
				row[i] = dup_node_map.get(int(row[i]), int(row[i]))

		wall_edges_new = rewrite_boundary_with_duplicate_nodes(wall_edges, y, y_te, y_tol, dup_node_map)
		farfield_edges_new = rewrite_boundary_with_duplicate_nodes(farfield_edges, y, y_te, y_tol, dup_node_map)

		wake_conn_upper = np.array([(a, b) for (a, b, _) in wake_upper], dtype=wall_conn_dtype)
		wake_conn_lower = np.array([(a, b) for (a, b, _) in wake_lower], dtype=wall_conn_dtype)

		replace_data_dataset(coords["CoordinateX"], x_new.astype(np.float64))
		replace_data_dataset(coords["CoordinateY"], y_new.astype(np.float64))
		replace_data_dataset(coords["CoordinateZ"], z_new.astype(np.float64))

		replace_data_dataset(zone[vol_name]["ElementConnectivity"], elem_conn_new.reshape(-1).astype(vol_conn_dtype, copy=False))
		n_volume = elem_conn_new.shape[0]
		zone[vol_name]["ElementRange"][" data"][...] = np.array([1, n_volume], dtype=vol_range_dtype)

		next_start = n_volume + 1
		next_start = update_section_connectivity(
			zone, farfield_name, farfield_edges_new.reshape(-1), next_start, conn_dtype=farfield_conn_dtype
		)
		next_start = update_section_connectivity(
			zone, wall_name, wall_edges_new.reshape(-1), next_start, conn_dtype=wall_conn_dtype
		)
		wake_ranges = None
		if write_wake_boundary:
			for name in ["wake", "wake_lower", "wake_upper"]:
				if name in zone and name not in [wall_name, farfield_name, vol_name]:
					del zone[name]

			for name in ["wake_lower", "wake_upper"]:
				if name not in zone:
					zone.copy(wall_name, name)
				set_cgns_name(zone[name], name)

			next_start = update_section_connectivity(
				zone, "wake_lower", wake_conn_lower.reshape(-1), next_start, conn_dtype=wall_conn_dtype
			)
			_ = update_section_connectivity(
				zone, "wake_upper", wake_conn_upper.reshape(-1), next_start, conn_dtype=wall_conn_dtype
			)
			wake_ranges = {
				"wake_lower": zone["wake_lower"]["ElementRange"][" data"][:],
				"wake_upper": zone["wake_upper"]["ElementRange"][" data"][:],
			}
		else:
			for name in ["wake", "wake_lower", "wake_upper"]:
				if name in zone:
					del zone[name]

		farfield_range = zone[farfield_name]["ElementRange"][" data"][:]
		wall_range = zone[wall_name]["ElementRange"][" data"][:]
		update_zonebc_tags(zone, farfield_range, wall_range, wake_ranges)
		ensure_base_wake_families(base, ["wake_lower", "wake_upper"], enabled=write_wake_boundary)

		zone_data = zone[" data"][:]
		zone_data[0, 0] = x_new.size
		zone[" data"][...] = zone_data

		if write_wake_boundary:
			write_wake_links_to_cgns(zone, wake_lower, wake_upper, wake_ranges, x_new, y_new, vol_range_dtype)
		elif "WakeLinkage" in zone:
			del zone["WakeLinkage"]

		print(f"Input mesh: {src}")
		print(f"Output mesh: {dst}")
		print(f"Volume section: {vol_name} ({n_per_elem} nodes/element)")
		print(f"Detected trailing edge: x = {x_te:.12g}, y = {y_te:.12g}")
		print(f"Wake interior faces split: {len(wake_edges)}")
		print(f"Wake nodes duplicated: {len(wake_nodes)}")
		print(f"Total nodes: {x.size} -> {x_new.size}")
		if write_wake_boundary:
			print(f"Wake lower edges: {len(wake_lower)}")
			print(f"Wake upper edges: {len(wake_upper)}")
			print("Wake linkage node: /Base/dom-1/WakeLinkage")
		else:
			print("Wake boundary entities: omitted")


def build_parser():
	parser = argparse.ArgumentParser(
		description=(
			"Split wake connectivity in a 2D CGNS mesh by duplicating centerline nodes "
			"from trailing edge to farfield so upper and lower cells are decoupled."
		)
	)
	parser.add_argument("input", help="Input CGNS mesh file")
	parser.add_argument(
		"-o",
		"--output",
		help="Output CGNS mesh file (default: <input>_wake.cgns)",
	)
	parser.add_argument(
		"--airfoil-file",
		default=None,
		help=(
			"Optional airfoil coordinate file (x y) used to determine trailing-edge location. "
			"If omitted, trailing edge is detected from wall boundary nodes."
		),
	)
	parser.add_argument(
		"--y-tol",
		type=float,
		default=1e-10,
		help="Tolerance for selecting wake-line nodes around y = y_TE",
	)
	parser.add_argument(
		"--x-tol",
		type=float,
		default=1e-12,
		help="Tolerance for selecting edges with x >= x_TE",
	)
	parser.add_argument(
		"--overwrite",
		action="store_true",
		help="Overwrite output file if it already exists",
	)
	parser.add_argument(
		"--omit-wake-boundary",
		action="store_true",
		help="Do not write wake section, ZoneBC/wake, or Base/wake family metadata",
	)
	return parser


def main():
	parser = build_parser()
	args = parser.parse_args()

	input_path = Path(args.input)
	if args.output:
		output_path = Path(args.output)
	else:
		output_path = input_path.with_name(f"{input_path.stem}_wake{input_path.suffix}")

	split_wake(
		mesh_in=str(input_path),
		mesh_out=str(output_path),
		airfoil_file=args.airfoil_file,
		y_tol=args.y_tol,
		x_tol=args.x_tol,
		overwrite=args.overwrite,
		write_wake_boundary=not args.omit_wake_boundary,
	)


if __name__ == "__main__":
	main()