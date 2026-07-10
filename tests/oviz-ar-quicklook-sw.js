const OVIZ_AR_QUICKLOOK_CACHE = "oviz-ar-quicklook-v2";

async function ovizArRangeResponse(request, response) {
  const range = request.headers.get("range");
  if (!range) {
    return response;
  }
  const match = /^bytes=(\d*)-(\d*)$/i.exec(range);
  if (!match) {
    return response;
  }
  const buffer = await response.arrayBuffer();
  const size = buffer.byteLength;
  let start = match[1] ? Number.parseInt(match[1], 10) : 0;
  let end = match[2] ? Number.parseInt(match[2], 10) : size - 1;
  if (!match[1] && match[2]) {
    const suffixLength = Math.max(0, Number.parseInt(match[2], 10) || 0);
    start = Math.max(size - suffixLength, 0);
    end = size - 1;
  }
  start = Math.max(0, Math.min(start, size - 1));
  end = Math.max(start, Math.min(end, size - 1));
  const headers = new Headers(response.headers);
  headers.set("Content-Range", `bytes ${start}-${end}/${size}`);
  headers.set("Accept-Ranges", "bytes");
  headers.set("Content-Length", String(end - start + 1));
  return new Response(buffer.slice(start, end + 1), {
    status: 206,
    statusText: "Partial Content",
    headers,
  });
}

self.addEventListener("install", (event) => {
  event.waitUntil(self.skipWaiting());
});

self.addEventListener("activate", (event) => {
  event.waitUntil(self.clients.claim());
});

self.addEventListener("fetch", (event) => {
  const url = new URL(event.request.url);
  if (!url.pathname.includes("/oviz-ar-quicklook/")) {
    return;
  }
  event.respondWith((async () => {
    const cache = await caches.open(OVIZ_AR_QUICKLOOK_CACHE);
    const cached = await cache.match(event.request, { ignoreSearch: false });
    if (cached) {
      return ovizArRangeResponse(event.request, cached);
    }
    return new Response("Missing Oviz AR snapshot.", {
      status: 404,
      headers: {
        "Content-Type": "text/plain; charset=utf-8",
        "Cache-Control": "no-store",
      },
    });
  })());
});
