// @ts-check
import { defineConfig } from 'astro/config';
import starlight from '@astrojs/starlight';
import starlightThemeNova from 'starlight-theme-nova';
import { readFileSync, existsSync } from 'node:fs';
import { fileURLToPath } from 'node:url';

const sidebarManifestPath = fileURLToPath(new URL('./build/sidebar.json', import.meta.url));
const referenceSidebar = existsSync(sidebarManifestPath)
	? JSON.parse(readFileSync(sidebarManifestPath, 'utf8'))
	: [];

// https://astro.build/config
export default defineConfig({
	integrations: [
		starlight({
			title: 'MIINT',
			description: 'Microbiome analytics for DuckDB.',
			plugins: [starlightThemeNova()],
			social: [
				{ icon: 'github', label: 'GitHub', href: 'https://github.com/the-miint/duckdb-miint' },
			],
			// "Edit this page on GitHub" link in the right column. Points at the
			// branch the docs are currently built from. Generated reference pages
			// (excluded by site/.gitignore) get an edit link too — the link 404s
			// for those, which is the right signal to "edit the C++ instead".
			editLink: {
				baseUrl: 'https://github.com/the-miint/duckdb-miint/edit/v1.5-variegata/site/',
			},
			// Show "Last updated <date>" in the page footer, sourced from the
			// file's git mtime.
			lastUpdated: true,
			// Show "next page" / "previous page" links in the page footer.
			pagination: true,
			sidebar: [
				{
					label: 'Getting started',
					items: [
						{ label: 'Installation', slug: 'getting-started/installation' },
						{ label: 'Updating', slug: 'getting-started/updating' },
					],
				},
				{
					label: 'Internals',
					autogenerate: { directory: 'internals' },
				},
				{
					label: 'Reference',
					items: referenceSidebar,
				},
			],
		}),
	],
});
