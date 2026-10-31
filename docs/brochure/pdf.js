// The brochure as a PDF (two A4 pages).
//
// Usage: node pdf.js brochure.html lumipore-brochure.pdf
//
// Requires Playwright (with Chromium).
const { chromium } = require('playwright');
const path = require('path');

(async () => {
  const [html, pdf] = process.argv.slice(2);
  const browser = await chromium.launch();
  const page = await browser.newPage();
  await page.goto('file://' + path.resolve(html));
  await page.waitForTimeout(500);
  await page.pdf({ path: pdf, format: 'A4', printBackground: true,
                   preferCSSPageSize: true });
  await browser.close();
})();
