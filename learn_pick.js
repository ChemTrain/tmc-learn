/* The Learn tab's next card, the same in every ChemTrain app (shared/frame.py puts it in every page).

   No rounds: one question after another until the set is mastered.
   - A set (a subcategory, or My cards): any card still to ask, at random; never one of the last few asked.
   - The All pool of an app: a card waits before it comes back, longer the better it is known
     (WAIT_ALL questions, by its level), and cards not yet mastered come up more often than mastered ones.
   Levels: three right answers master a card (Not studied → Learning → Familiar → Mastered);
   a wrong answer moves it back one step, never below Learning.
   The other quiz tabs (Most acidic, Which is which, Removed?, Plan, the Reactivity quizzes, Count, Name) work the
   same way through deal/step/mark: no rounds of eight, a question counts as done after three right answers. */
window.CTLearn = (function () {
  var WAIT_SET = 3;                    // a set: not one of the last 3 cards asked (fewer in a tiny set)
  var WAIT_ALL = [10, 10, 30, 80];     // All: questions before a card can come back, by its level
  var MASTERED_WEIGHT = 0.25;          // All, Ask "All cards": a mastered card is drawn a quarter as often

  function isAll(set) { return !!set && (set.id === 'all' || /^all-/.test(String(set.id))); }
  function after(before, ok) { return ok ? Math.min(3, before + 1) : Math.max(1, before - 1); }

  // o: { cards, order, lvl, ask ('all' | anything else: not mastered only), shuffle, all, round, not, key }
  // round.seen maps a card's key (its id) to the question number it was last asked at. Returns the card, or null when done.
  function next(o) {
    var key = o.key || function (c) { return c.id; };
    var lvl = function (c) { return o.lvl(key(c)); }, round = o.round, seen = round.seen || (round.seen = new Map());
    var pool = o.cards.filter(function (c) { return o.ask === 'all' || lvl(c) < 3; });
    if (!pool.length) return null;
    if (o.not && pool.length > 1) pool = pool.filter(function (c) { return c !== o.not; });
    var n = round.queue.length, pick;
    if (!o.shuffle) {  // Set order: the next card after the last one asked
      var order = o.order.filter(function (c) { return pool.indexOf(c) >= 0; });
      var last = -1;
      order.forEach(function (c, i) { if (seen.has(key(c)) && (last < 0 || seen.get(key(c)) > seen.get(key(order[last])))) last = i; });
      pick = order[(last + 1) % order.length];
    } else {
      var wait = function (c) { return Math.min(pool.length - 1, o.all ? WAIT_ALL[lvl(c)] : WAIT_SET); };
      var since = function (c) { return seen.has(key(c)) ? n - seen.get(key(c)) : Infinity; };
      var free = pool.filter(function (c) { return since(c) > wait(c); });
      if (!free.length) {  // every card asked too recently: the one closest to its turn
        var best = Infinity;
        pool.forEach(function (c) { var d = wait(c) - since(c); if (d < best) best = d; });
        free = pool.filter(function (c) { return wait(c) - since(c) === best; });
      }
      var w = function (c) { return o.all && lvl(c) === 3 ? MASTERED_WEIGHT : 1; };
      var total = free.reduce(function (s, c) { return s + w(c); }, 0), r = Math.random() * total;
      pick = free[free.length - 1];
      for (var i = 0; i < free.length; i++) { r -= w(free[i]); if (r < 0) { pick = free[i]; break; } }
    }
    seen.set(key(pick), n);
    return pick;
  }

  // A quiz tab as a stream: st.deck holds the questions asked so far, st.i the one on screen.
  // q: { pool, key, rec (key -> right answers so far; 3 = done), ask, all, not }. Each returns the question, or null.
  function pickInto(st, q, not) {
    var c = next({ cards: q.pool, key: q.key, lvl: function (k) { return Math.min(3, +q.rec[k] || 0); }, ask: q.ask,
      shuffle: true, all: q.all, round: { queue: st.deck, seen: st.seenAt }, not: not });
    if (c) st.deck[st.i] = c;
    return c;
  }
  function deal(st, q, not) { st.deck = []; st.i = 0; st.results = []; st.retry = false; st.seenAt = new Map(); return pickInto(st, q, not); }
  function step(st, q) { if (!st.seenAt) return deal(st, q); st.i++; return pickInto(st, q); }
  function swap(st, q) { return st.seenAt ? pickInto(st, q, st.deck[st.i]) : deal(st, q); }  // another question in place of the one on screen
  function mark(rec, k, ok) { rec[k] = after(+rec[k] || 0, ok); }
  var done = function (rec, k) { return (+rec[k] || 0) >= 3; };

  return { next: next, isAll: isAll, after: after, deal: deal, step: step, swap: swap, mark: mark, done: done, WAIT_SET: WAIT_SET, WAIT_ALL: WAIT_ALL };
})();
