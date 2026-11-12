/* How a ChemTrain account's saved data merges across devices. Shared by the pages (inlined by shared/frame.py,
   which drops the export line) and the server (net/src/worker.js imports it).

   A synced document maps a localStorage key to { v: its value, t: when that device last changed it }. */

// What to do with a key: skip it, keep the newer value (lww), keep the higher number per entry (max), keep every
// entry (union), or keep the play with more answers (progress).
function KIND(k) {
  if (!/^(pka2|rxn|pg|tmc1|rl|rad|pr|ctd)\./.test(k)) return 'skip';                    // not an app's record
  if (/\.round\.|\.pt$|^ctd\.past\.|^ctd\.set$/.test(k)) return 'skip';                 // a round or set under way, an open panel
  if (/^ctd\.daily\./.test(k)) return 'progress';                                        // a Daily quiz: its questions and answers
  if (/\.levels(\.|$)|\.seen$|^(ctd\.history|ctd\.late|ctd\.beast|pg\.plan|pg\.removed|rxn\.which|pka2\.acidic|tmc1\.named)$/.test(k)) return 'max';
  if (/\.stars$/.test(k)) return 'union';                                                // starred cards
  return 'lww';                                                                          // settings, picks, custom cards
}

const isObj = v => v && typeof v === 'object' && !Array.isArray(v);

// one key: the record kept when `a` (what is stored) meets `b` (what arrives)
function mergeRec(kind, a, b) {
  if (!a) return b;
  if (!b) return a;
  const newer = (b.t || 0) >= (a.t || 0) ? b : a, older = newer === b ? a : b, t = Math.max(a.t || 0, b.t || 0);
  if (kind === 'max' && isObj(a.v) && isObj(b.v)) {
    const v = { ...older.v, ...newer.v };
    for (const k of Object.keys(v)) {
      const x = a.v[k], y = b.v[k];
      if (typeof x === 'number' && typeof y === 'number') v[k] = Math.max(x, y);
    }
    return { v, t };
  }
  if (kind === 'union') {
    if (Array.isArray(a.v) && Array.isArray(b.v)) return { v: [...new Set([...a.v, ...b.v])], t };
    if (isObj(a.v) && isObj(b.v)) return { v: { ...older.v, ...newer.v }, t };
  }
  if (kind === 'progress' && isObj(a.v) && isObj(b.v)) {
    const na = Array.isArray(a.v.answers) ? a.v.answers.length : 0, nb = Array.isArray(b.v.answers) ? b.v.answers.length : 0;
    return na === nb ? { v: newer.v, t } : na > nb ? { v: a.v, t } : { v: b.v, t };
  }
  return { v: newer.v, t };
}

// the stored document after `changes` from a device
function mergeDocs(doc, changes) {
  const out = { ...doc };
  for (const [k, rec] of Object.entries(changes || {})) {
    const kind = KIND(k);
    if (kind === 'skip' || !rec || !('v' in rec)) continue;
    out[k] = mergeRec(kind, out[k], { v: rec.v, t: Number(rec.t) || 0 });
  }
  return out;
}

export { KIND, mergeDocs };
