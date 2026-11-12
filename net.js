/* ChemTrain's account on this device, and its progress kept in step across devices. Inlined into every ChemTrain
   page by shared/frame.py, right after sync_rules.js (which gives KIND), before the page's own script.

   An account is { id, key, pseudo } under localStorage 'ct.account'. Without one, nothing here does anything.
   With one, the page's saved records (see KIND) go to the server and the other devices' changes come back:
   at load, when the page is hidden, when the connection returns, and every minute while something changed.
   The server merges (shared/sync_rules.js) and answers with the whole document, which is written back here.
   A change that arrives after the page has started, before anything was touched, reloads it once. */
(function () {
  'use strict';
  const DEFAULT_API = 'https://api.chemtrain.workers.dev';
  const json = k => { try { const v = localStorage.getItem(k); return v ? JSON.parse(v) : null; } catch (e) { return null; } };
  const save = (k, v) => { try { v === null ? localStorage.removeItem(k) : localStorage.setItem(k, JSON.stringify(v)); } catch (e) {} };
  const API = json('ct.api') || DEFAULT_API;
  const account = () => json('ct.account');

  async function api(path, body, opts) {
    const init = body ? { method: 'POST', headers: { 'content-type': 'application/json' }, body: JSON.stringify(body), ...opts } : { method: 'GET' };
    const r = await fetch(API + path, init);
    const d = await r.json().catch(() => ({}));
    if (!r.ok) throw Object.assign(new Error(d.error || 'HTTP ' + r.status), { status: r.status, code: d.error || 'http_' + r.status });
    return d;
  }
  const authed = (path, body, opts) => { const a = account(); return a ? api(path, { key: a.key, ...body }, opts) : Promise.reject(Object.assign(new Error('no account'), { code: 'no_account' })); };

  /* ---------- progress across devices ---------- */
  const hash = s => { let h = 5381; for (let i = 0; i < s.length; i++) h = (h * 33 + s.charCodeAt(i)) | 0; return h.toString(36); };
  function keysToSync() {
    const cut = new Date(Date.now() - 60 * 864e5).toISOString().slice(0, 10), out = [];  // Daily quizzes: the last 60 days
    for (let i = 0; i < localStorage.length; i++) {
      const k = localStorage.key(i), kind = KIND(k);
      if (kind === 'skip' || (kind === 'progress' && k.slice('ctd.daily.'.length) < cut)) continue;
      out.push(k);
    }
    return out;
  }
  function localChanges() {
    const snap = json('ct.sync.snap') || {}, changes = {}, now = Date.now();
    for (const k of keysToSync()) {
      const raw = localStorage.getItem(k);
      if (raw == null || snap[k] === hash(raw)) continue;
      try { changes[k] = { v: JSON.parse(raw), t: now }; } catch (e) {}
    }
    return changes;
  }
  let booted = document.readyState !== 'loading', touched = false, busy = null;
  document.addEventListener('DOMContentLoaded', () => { booted = true; });
  for (const ev of ['pointerdown', 'keydown']) addEventListener(ev, () => { touched = true; }, { capture: true, once: true });

  function sync(opts) {
    const a = account();
    if (!a || busy) return busy || Promise.resolve([]);
    busy = (async () => {
      const changes = localChanges(), body = { changes, have: json('ct.sync.updated') || 0 };
      const text = JSON.stringify(body), init = opts && opts.leaving && text.length < 60000 ? { keepalive: true } : undefined;
      const res = await authed('/sync', body, init);
      const snap = json('ct.sync.snap') || {}, changed = [];
      if (res.doc) {
        for (const [k, rec] of Object.entries(res.doc)) {
          if (KIND(k) === 'skip') continue;
          const raw = JSON.stringify(rec.v);
          if (localStorage.getItem(k) !== raw) { localStorage.setItem(k, raw); changed.push(k); }
          snap[k] = hash(raw);
        }
        save('ct.sync.updated', res.updated);
      }
      save('ct.sync.snap', snap); save('ct.sync.at', Date.now());
      if (changed.length) {
        dispatchEvent(new CustomEvent('ct-sync', { detail: { changed } }));
        if (booted && !touched && !document.hidden && !(opts && opts.noReload)) location.reload();
      }
      return changed;
    })().catch(() => []).finally(() => { busy = null; });
    return busy;
  }
  const dirty = () => Object.keys(localChanges()).length > 0;

  /* ---------- Daily scores ---------- */
  // today's score goes up when there is an account; otherwise it waits for one
  async function submitScore(day, score, total) {
    const pend = json('ct.score.pending') || {}; pend[day] = { score, total }; save('ct.score.pending', pend);
    return account() ? flushScores() : null;
  }
  async function flushScores() {
    if (!account()) return null;
    const pend = json('ct.score.pending') || {}; let last = null;
    for (const day of Object.keys(pend)) {
      try { last = await authed('/score', { day, ...pend[day] }); delete pend[day]; save('ct.score.' + day, last); }
      catch (e) { if (e.status === 400) delete pend[day]; }
    }
    save('ct.score.pending', pend);
    return last;
  }
  const scoreOf = day => json('ct.score.' + day);

  /* ---------- the account ---------- */
  function setAccount(a) { save('ct.account', a); save('ct.sync.snap', null); save('ct.sync.updated', null); }
  async function join(pseudo) { const r = await api('/join', { pseudo }); setAccount(r); await Promise.all([sync({ noReload: true }), flushScores()]); return r; }
  async function login(key) { const r = await api('/login', { key }); setAccount({ id: r.id, key: String(key).toUpperCase().replace(/[^0-9A-Z]/g, '').replace(/O/g, '0').replace(/[IL]/g, '1').match(/.{1,4}/g).join('-'), pseudo: r.pseudo }); await Promise.all([sync({ noReload: true }), flushScores()]); return r; }
  async function rename(pseudo) { const r = await authed('/pseudo', { pseudo }); save('ct.account', { ...account(), pseudo: r.pseudo }); return r; }
  function logout() { for (const k of ['ct.account', 'ct.sync.snap', 'ct.sync.updated', 'ct.sync.at']) save(k, null); }
  async function remove() { await authed('/delete', {}); logout(); }

  addEventListener('visibilitychange', () => { if (document.hidden) sync({ leaving: true }); });
  addEventListener('online', () => { sync(); flushScores(); });
  setInterval(() => { if (account() && !document.hidden && dirty()) sync(); }, 60000);

  window.CT = { API, api, authed, account, join, login, logout, rename, remove, sync, submitScore, flushScores, scoreOf, KIND, hash };
  if (account()) { sync(); if (Object.keys(json('ct.score.pending') || {}).length) flushScores(); }
})();
