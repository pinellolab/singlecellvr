/*global THREE, AFRAME, Utils*/

AFRAME.registerComponent('mimic-rotation', {

    schema: {
        masterId: {type: "string"}
    },

    tick: function() {
        const masterEl = document.getElementById(this.data.masterId);
        const [xR, yR, zR] = Object.values(masterEl.object3D.rotation);
        this.el.object3D.rotation.set(
            xR, yR, zR
        );
    }
});
