---
layout: page
title: About
permalink: /about/
---


## About

This is a software developed by the [Steven Kelly Lab](http://www.stevekellylab.com/)

You can find the source code for OrthoFinder at GitHub:
[OrthoFinder](https://github.com/OrthoFinder/OrthoFinder)

### Why OrthoFinder?





### Meet&nbsp;the&nbsp;Team {#team}

#### Current&nbsp;Team&nbsp;Member {#team-current}

<div class="team-grid">
{% assign current_team = site.data.team | where: "status", "Current" %}
{% for member in current_team %}
  {% include team-card.html member=member %}
{% endfor %}
</div>

---

#### Former&nbsp;Team&nbsp;Member {#team-former}

<div class="team-grid">
{% assign former_team = site.data.team | where: "status", "Former" %}
{% for member in former_team %}
  {% include team-card.html member=member %}
{% endfor %}
</div>
