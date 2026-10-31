#define PY_SSIZE_T_CLEAN
#include <Python.h>

#include <algorithm>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "readset.h"

struct PythonError {};

static void check(int status) {
	if (status < 0) {
		throw PythonError();
	}
}

static PyObject* check(PyObject* object) {
	if (object == nullptr) {
		throw PythonError();
	}
	return object;
}

struct Score {
	int first;
	int second;
	int quality;
};

static bool lower_score(const Score& left, const Score& right) {
	if (left.first != right.first) return left.first < right.first;
	if (left.second != right.second) return left.second < right.second;
	return left.quality < right.quality;
}

struct QueueEntry {
	Score score;
	int item;
};

class PriorityQueue {
public:
	std::vector<QueueEntry> heap;
	std::vector<int> positions;

	explicit PriorityQueue(size_t read_count) : positions(read_count, -1) {}

	void swap(size_t first, size_t second) {
		std::swap(heap[first], heap[second]);
		positions[heap[first].item] = first;
		positions[heap[second].item] = second;
	}

	bool lower(size_t first, size_t second) const {
		return lower_score(heap[first].score, heap[second].score);
	}

	void sift_up(size_t index) {
		while (index > 0) {
			size_t parent = (index - 1) / 2;
			if (!lower(parent, index)) break;
			swap(parent, index);
			index = parent;
		}
	}

	void sift_down(size_t index) {
		while (true) {
			size_t left = index * 2 + 1;
			size_t right = left + 1;
			if (right < heap.size()) {
				size_t child = lower(left, right) ? right : left;
				if (!lower(index, child)) break;
				swap(index, child);
				index = child;
			} else if (left < heap.size()) {
				if (!lower(index, left)) break;
				swap(index, left);
				index = left;
			} else {
				break;
			}
		}
	}

	void push(QueueEntry entry) {
		size_t index = heap.size();
		heap.push_back(entry);
		positions[entry.item] = index;
		sift_up(index);
	}

	QueueEntry pop() {
		QueueEntry first = heap.front();
		QueueEntry last = heap.back();
		heap.pop_back();
		positions[first.item] = -1;
		if (!heap.empty()) {
			heap[0] = last;
			positions[last.item] = 0;
			sift_down(0);
		}
		return first;
	}

	void change_score(int item, Score score) {
		size_t index = positions[item];
		Score old = heap[index].score;
		heap[index].score = score;
		if (lower_score(old, score)) {
			sift_up(index);
		} else {
			sift_down(index);
		}
	}
};

struct NativeIndex {
	std::unordered_map<int, long> positions;
	std::vector<PyObject*> readers;
	std::vector<Score> scores;
};

static long index_of(const NativeIndex& indices, int position) {
	auto entry = indices.positions.find(position);
	return entry == indices.positions.end() ? -1 : entry->second;
}

struct Coverage {
	std::vector<int> values;
	explicit Coverage(size_t size) : values(size, 0) {}
	long maximum(int begin, int end) const {
		if (begin < 0 || begin >= end || static_cast<size_t>(end) > values.size()) {
			throw std::out_of_range("coverage range " + std::to_string(begin) + " " + std::to_string(end) + " " + std::to_string(values.size()));
		}
		return *std::max_element(values.begin() + begin, values.begin() + end);
	}
	void add(int begin, int end) {
		for (int index = begin; index < end; ++index) {
			++values[index];
		}
	}
};

struct Components {
	std::vector<int> parents;
	explicit Components(size_t size) {
		parents.reserve(size);
		for (size_t index = 0; index < size; ++index) parents.push_back(index);
	}
	int find(int index) {
		if (index < 0 || static_cast<size_t>(index) >= parents.size()) {
			throw std::out_of_range("component index " + std::to_string(index));
		}
		int root = index;
		while (parents[root] != root) root = parents[root];
		while (parents[index] != index) {
			int next = parents[index];
			parents[index] = root;
			index = next;
		}
		return root;
	}
	void merge(int first, int second) {
		first = find(first);
		second = find(second);
		if (first != second) {
			int small = std::min(first, second);
			int large = std::max(first, second);
			parents[large] = small;
		}
	}
};

static Score read_score(ReadSet* reads, int index, const NativeIndex& indices) {
	Read* read = reads->get(index);
	int minimum_quality = -1;
	int good_score = 0;
	long first = -1;
	long last = -1;
	for (int i = 0; i < read->getVariantCount(); ++i) {
		int quality = read->getVariantQuality(i);
		if (i == 0) minimum_quality = quality;
		else minimum_quality = std::min(minimum_quality, quality);
		long position = index_of(indices, read->getPosition(i));
		if (position >= 0) {
			if (good_score == 0) first = position;
			last = position;
			++good_score;
		}
	}
	if (good_score == 0) {
		throw std::out_of_range("no positions for read " + std::to_string(index));
	}
	long span = last - first + 1;
	int bad_score = span == good_score ? 0 : static_cast<int>(span - good_score);
	return {good_score - bad_score, good_score - bad_score, minimum_quality};
}

static PriorityQueue make_queue(ReadSet* reads, const NativeIndex& indices, PyObject* read_indices) {
	PriorityQueue queue(reads->size());
	PyObject* iterator = check(PyObject_GetIter(read_indices));
	PyObject* item;
	while ((item = PyIter_Next(iterator)) != nullptr) {
		long index = PyLong_AsLong(item);
		Py_DECREF(item);
		if (index < 0 && PyErr_Occurred()) throw PythonError();
		queue.push({indices.scores[index], static_cast<int>(index)});
	}
	if (PyErr_Occurred()) throw PythonError();
	Py_DECREF(iterator);
	return queue;
}

static long maximum_coverage(const Coverage& coverages, int begin, int end) {
	return coverages.maximum(begin, end);
}

static void add_coverage(Coverage& coverages, int begin, int end) {
	coverages.add(begin, end);
}

static void update_set(PyObject* set, PyObject* iterable) {
	PyObject* result = check(PyObject_CallMethod(set, "update", "O", iterable));
	Py_DECREF(result);
}

static void add_index(PyObject* set, int index) {
	PyObject* value = check(PyLong_FromLong(index));
	check(PySet_Add(set, value));
	Py_DECREF(value);
}

static bool contains_index(PyObject* set, int index) {
	PyObject* value = check(PyLong_FromLong(index));
	int present = PySet_Contains(set, value);
	Py_DECREF(value);
	if (present < 0) throw PythonError();
	return present != 0;
}

static void discard_index(PyObject* set, int index) {
	PyObject* value = check(PyLong_FromLong(index));
	check(PySet_Discard(set, value));
	Py_DECREF(value);
}

static void merge_position(Components& finder, const NativeIndex& indices, int first, int second) {
	finder.merge(index_of(indices, first), index_of(indices, second));
}

static int find_position(Components& finder, int position) {
	return finder.find(position);
}

static void merge_read(Components& finder, const NativeIndex& indices, Read* read) {
	for (int i = 1; i < read->getVariantCount(); ++i) {
		merge_position(finder, indices, read->getPosition(0), read->getPosition(i));
	}
}

struct Slice {
	PyObject* selected;
	PyObject* violating;
};

static Slice select_slice(ReadSet* reads, Coverage& coverages, int max_coverage,
	const NativeIndex& indices, PyObject* undecided) {
	PriorityQueue queue = make_queue(reads, indices, undecided);
	PyObject* selected = check(PySet_New(nullptr));
	PyObject* violating = check(PySet_New(nullptr));
	std::unordered_set<int> covered;
	while (!queue.heap.empty()) {
		QueueEntry entry = queue.pop();
		Read* read = reads->get(entry.item);
		bool covers_new = false;
		std::unordered_set<int> new_positions;
		for (int i = 0; i < read->getVariantCount(); ++i) {
			int position = read->getPosition(i);
			if (covered.count(position) != 0) {
				continue;
			}
			covers_new = true;
			new_positions.insert(position);
		}
		int begin = index_of(indices, read->getPosition(0));
		int end = index_of(indices, read->getPosition(read->getVariantCount() - 1)) + 1;
		if (maximum_coverage(coverages, begin, end) >= max_coverage) {
			add_index(violating, entry.item);
		} else if (covers_new) {
			add_coverage(coverages, begin, end);
			add_index(selected, entry.item);
			PyObject* updates = check(PySet_New(nullptr));
			for (int position : new_positions) {
				covered.insert(position);
				update_set(updates, indices.readers[index_of(indices, position)]);
			}
			PyObject* difference = check(PyObject_CallMethod(updates, "difference", "O", selected));
			PyObject* iterator = check(PyObject_GetIter(difference));
			PyObject* item;
			while ((item = PyIter_Next(iterator)) != nullptr) {
				int index = PyLong_AsLong(item);
				Py_DECREF(item);
				if (queue.positions[index] >= 0) {
					Score old = queue.heap[queue.positions[index]].score;
					Read* affected = reads->get(index);
					int uncovered = 0;
					for (int i = 0; i < affected->getVariantCount(); ++i) {
						uncovered += new_positions.count(affected->getPosition(i)) == 0;
					}
					queue.change_score(index, {old.first - uncovered, old.second, old.quality});
				}
			}
			if (PyErr_Occurred()) throw PythonError();
			Py_DECREF(iterator);
			Py_DECREF(difference);
			Py_DECREF(updates);
		}
	}
	return {selected, violating};
}

static PyObject* select_reads(ReadSet* reads, int max_coverage, PyObject* positions,
	PyObject* indices, PyObject* variant_to_reads) {
	NativeIndex native;
	Py_ssize_t offset = 0;
	PyObject* key;
	PyObject* value;
	while (PyDict_Next(indices, &offset, &key, &value)) {
		long position = PyLong_AsLong(key);
		long index = PyLong_AsLong(value);
		native.positions[position] = index;
		if (static_cast<size_t>(index) >= native.readers.size()) {
			native.readers.resize(index + 1, nullptr);
		}
		native.readers[index] = check(PyDict_GetItemWithError(variant_to_reads, value));
	}
	native.scores.reserve(reads->size());
	for (size_t index = 0; index < reads->size(); ++index) {
		native.scores.push_back(read_score(reads, index, native));
	}
	Coverage coverages(PyObject_Length(positions));
	PyObject* selected = check(PySet_New(nullptr));
	PyObject* range = check(PyObject_CallFunction((PyObject*)&PyRange_Type, "i", static_cast<int>(reads->size())));
	PyObject* undecided = check(PySet_New(range));
	Py_DECREF(range);
	while (PySet_Size(undecided) > 0) {
		Slice slice = select_slice(reads, coverages, max_coverage, native, undecided);
		update_set(selected, slice.selected);
		PyObject* next = check(PyNumber_InPlaceSubtract(undecided, slice.selected));
		Py_DECREF(undecided);
		undecided = next;
		next = check(PyNumber_InPlaceSubtract(undecided, slice.violating));
		Py_DECREF(undecided);
		undecided = next;
		Components finder(native.positions.size());
		PyObject* iterator = check(PyObject_GetIter(slice.selected));
		PyObject* item;
		while ((item = PyIter_Next(iterator)) != nullptr) {
			int index = PyLong_AsLong(item);
			Py_DECREF(item);
			merge_read(finder, native, reads->get(index));
		}
		if (PyErr_Occurred()) throw PythonError();
		Py_DECREF(iterator);
		PyObject* bridging = check(PySet_New(nullptr));
		PriorityQueue queue = make_queue(reads, native, undecided);
		while (!queue.heap.empty()) {
			QueueEntry entry = queue.pop();
			Read* read = reads->get(entry.item);
			std::unordered_set<int> blocks;
			for (int i = 0; i < read->getVariantCount(); ++i) {
				blocks.insert(find_position(finder, index_of(native, read->getPosition(i))));
			}
			int begin = index_of(native, read->getPosition(0));
			int end = index_of(native, read->getPosition(read->getVariantCount() - 1)) + 1;
			if (maximum_coverage(coverages, begin, end) >= max_coverage) {
				discard_index(undecided, entry.item);
				continue;
			}
			if (blocks.size() < 2) continue;
			add_index(bridging, entry.item);
			add_index(selected, entry.item);
			add_coverage(coverages, begin, end);
			discard_index(undecided, entry.item);
			merge_read(finder, native, read);
		}
		Py_DECREF(bridging);
		Py_DECREF(slice.selected);
		Py_DECREF(slice.violating);
	}
	Py_DECREF(undecided);
	return selected;
}

static PyObject* select(PyObject*, PyObject* args) {
	try {
		PyObject* readset_object;
		PyObject* positions;
		PyObject* indices;
		PyObject* variant_to_reads;
		int max_coverage;
		if (!PyArg_ParseTuple(args, "OiOOO", &readset_object, &max_coverage, &positions, &indices, &variant_to_reads)) {
			return nullptr;
		}
		PyObject* core = check(PyImport_ImportModule("whatshap.core"));
		PyObject* type = check(PyObject_GetAttrString(core, "ReadSet"));
		if (!PyObject_TypeCheck(readset_object, reinterpret_cast<PyTypeObject*>(type))) {
			PyErr_SetString(PyExc_TypeError, "expected whatshap.core.ReadSet");
			return nullptr;
		}
		struct PythonReadSet {
			PyObject_HEAD
			ReadSet* thisptr;
		};
		ReadSet* reads = reinterpret_cast<PythonReadSet*>(readset_object)->thisptr;
		return select_reads(reads, max_coverage, positions, indices, variant_to_reads);
	} catch (PythonError) {
		return nullptr;
	} catch (const std::exception& error) {
		PyErr_SetString(PyExc_RuntimeError, error.what());
		return nullptr;
	}
}

struct PythonReadSetView {
	PyObject_HEAD
	ReadSet* reads;
};

static ReadSet* get_readset(PyObject* object) {
	PyObject* core = check(PyImport_ImportModule("whatshap.core"));
	PyObject* type = check(PyObject_GetAttrString(core, "ReadSet"));
	if (!PyObject_TypeCheck(object, reinterpret_cast<PyTypeObject*>(type))) {
		PyErr_SetString(PyExc_TypeError, "expected whatshap.core.ReadSet");
		throw PythonError();
	}
	return reinterpret_cast<PythonReadSetView*>(object)->reads;
}

static void write_u32(std::string& output, uint32_t value) {
	for (int index = 0; index < 4; ++index) {
		output.push_back(static_cast<char>((value >> (index * 8)) & 0xff));
	}
}

static void write_int(std::string& output, int value) {
	write_u32(output, static_cast<uint32_t>(value));
}

static void write_text(std::string& output, const std::string& value) {
	write_u32(output, value.size());
	output.append(value);
}

struct BufferReader {
	const char* data;
	size_t size;
	size_t offset;
	uint32_t u32() {
		if (offset + 4 > size) throw std::runtime_error("truncated read-set data");
		uint32_t value = 0;
		for (int index = 0; index < 4; ++index) {
			value |= static_cast<uint32_t>(static_cast<unsigned char>(data[offset++])) << (index * 8);
		}
		return value;
	}
	int integer() {
		return static_cast<int32_t>(u32());
	}
	std::string text() {
		uint32_t length = u32();
		if (offset + length > size) throw std::runtime_error("truncated read-set text");
		std::string result(data + offset, length);
		offset += length;
		return result;
	}
};

static PyObject* serialize(PyObject*, PyObject* args) {
	try {
		PyObject* object;
		if (!PyArg_ParseTuple(args, "O", &object)) return nullptr;
		ReadSet* reads = get_readset(object);
		std::string output("C3R1", 4);
		uint32_t count = 0;
		for (size_t index = 0; index < reads->size(); ++index) {
			count += reads->get(index)->getVariantCount() >= 2;
		}
		write_u32(output, count);
		for (size_t index = 0; index < reads->size(); ++index) {
			Read* read = reads->get(index);
			if (read->getVariantCount() < 2) continue;
			read->sortVariants();
			write_text(output, read->getName());
			write_u32(output, read->getMapqs().size());
			for (int mapq : read->getMapqs()) write_int(output, mapq);
			write_int(output, read->getSourceID());
			write_int(output, read->getSampleID());
			write_int(output, read->getReferenceStart());
			write_int(output, read->getReferenceEnd());
			write_int(output, read->getHPTag());
			write_int(output, read->getPSTag());
			write_int(output, read->isSupplementary());
			write_int(output, read->isReverse());
			write_text(output, read->getBXTag());
			write_text(output, read->getChromosome());
			write_text(output, read->getReadSubAlignmentId());
			write_u32(output, read->getVariantCount());
			for (int variant = 0; variant < read->getVariantCount(); ++variant) {
				write_int(output, read->getPosition(variant));
				write_int(output, read->getAllele(variant));
				write_int(output, read->getVariantQuality(variant));
			}
		}
		return PyBytes_FromStringAndSize(output.data(), output.size());
	} catch (PythonError) {
		return nullptr;
	} catch (const std::exception& error) {
		PyErr_SetString(PyExc_RuntimeError, error.what());
		return nullptr;
	}
}

static PyObject* deserialize(PyObject*, PyObject* args) {
	try {
		PyObject* blobs;
		if (!PyArg_ParseTuple(args, "O", &blobs)) return nullptr;
		PyObject* core = check(PyImport_ImportModule("whatshap.core"));
		PyObject* type = check(PyObject_GetAttrString(core, "ReadSet"));
		PyObject* object = check(PyObject_CallFunctionObjArgs(type, nullptr));
		ReadSet* reads = get_readset(object);
		PyObject* items = check(PySequence_Fast(blobs, "expected serialized read sets"));
		Py_ssize_t count = PySequence_Fast_GET_SIZE(items);
		for (Py_ssize_t index = 0; index < count; ++index) {
			PyObject* blob = PySequence_Fast_GET_ITEM(items, index);
			if (!PyBytes_Check(blob)) {
				PyErr_SetString(PyExc_TypeError, "expected bytes containing a read set");
				return nullptr;
			}
			BufferReader reader{PyBytes_AS_STRING(blob), static_cast<size_t>(PyBytes_GET_SIZE(blob)), 0};
			if (reader.size < 8 || std::string(reader.data, 4) != "C3R1") {
				throw std::runtime_error("invalid read-set serialization");
			}
			reader.offset = 4;
			uint32_t read_count = reader.u32();
			for (uint32_t read_index = 0; read_index < read_count; ++read_index) {
				std::string name = reader.text();
				uint32_t mapq_count = reader.u32();
				std::vector<int> mapqs;
				mapqs.reserve(mapq_count);
				for (uint32_t mapq = 0; mapq < mapq_count; ++mapq) mapqs.push_back(reader.integer());
				int source = reader.integer();
				int sample = reader.integer();
				int reference_start = reader.integer();
				int reference_end = reader.integer();
				int hp = reader.integer();
				int ps = reader.integer();
				bool supplementary = reader.integer() != 0;
				bool reverse = reader.integer() != 0;
				std::string bx = reader.text();
				std::string chromosome = reader.text();
				std::string subalignment = reader.text();
				Read* read = new Read(name, mapqs.empty() ? 0 : mapqs[0], source, sample,
					reference_start, bx, hp, ps, chromosome, subalignment,
					supplementary, reference_end, reverse);
				for (size_t mapq = 1; mapq < mapqs.size(); ++mapq) read->addMapq(mapqs[mapq]);
				uint32_t variant_count = reader.u32();
				for (uint32_t variant = 0; variant < variant_count; ++variant) {
					read->addVariant(reader.integer(), reader.integer(), reader.integer());
				}
				reads->add(read);
			}
		}
		Py_DECREF(items);
		reads->sort();
		return object;
	} catch (PythonError) {
		return nullptr;
	} catch (const std::exception& error) {
		PyErr_SetString(PyExc_RuntimeError, error.what());
		return nullptr;
	}
}

static PyMethodDef methods[] = {
	{"select", select, METH_VARARGS, "Select phase-informative reads"},
	{"serialize", serialize, METH_VARARGS, "Serialize a read set"},
	{"deserialize", deserialize, METH_VARARGS, "Deserialize read sets"},
	{nullptr, nullptr, 0, nullptr},
};

static struct PyModuleDef module = {
	PyModuleDef_HEAD_INIT,
	"clair3_readselect",
	nullptr,
	-1,
	methods,
};

PyMODINIT_FUNC PyInit_clair3_readselect() {
	return PyModule_Create(&module);
}
